#!/usr/bin/env python3
"""
Batch donor_atoms re-enrichment for BiometalDB using multiprocessing.
Re-processes ALL complexes with the FIXED ligand-splitting logic (GetMolFrags).
Designed for cloud instances (DO droplets) with 4+ vCPUs.

Usage:
    python batch_enrich_do.py [--db PATH] [--workers N] [--batch-size N] [--dry-run]
"""
import sqlite3
import json
import os
import sys
import argparse
import time
import signal
from multiprocessing import Pool, cpu_count
from collections import Counter

from rdkit import Chem

# ── Donor identification logic (from enrich_donor_atoms.py, fixed) ──

DONOR_ELEMENTS = {"N", "O", "S", "P", "Cl", "Br", "I", "F", "Se", "As"}

COORD_NUM = {
    ("Ru", 2): (4, 6), ("Ru", 3): (6, 6),
    ("Ir", 3): (6, 6),
    ("Rh", 1): (4, 5), ("Rh", 3): (6, 6),
    ("Os", 2): (6, 6), ("Os", 3): (6, 6), ("Os", 4): (6, 6),
    ("Re", 1): (6, 6), ("Re", 2): (6, 6), ("Re", 3): (6, 6), ("Re", 5): (6, 7),
}


def safe_mol_from_smiles(smiles):
    """Parse SMILES with maximum safety - no sanitization, catch all errors."""
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        return mol
    except Exception:
        return None


def split_ligands(smiles):
    """Split SMILES into individual ligand fragments. Skip single-atom counterions."""
    mol = safe_mol_from_smiles(smiles)
    if mol is None:
        return []
    try:
        frags = Chem.GetMolFrags(mol, asMols=True)
    except Exception:
        return []
    ligands = []
    for frag in frags:
        try:
            atoms = [a for a in frag.GetAtoms()]
            symbols = [a.GetSymbol() for a in atoms]
            if len(atoms) == 1:
                continue  # skip all single-atom fragments
            lig_smi = Chem.MolToSmiles(frag)
            if lig_smi:
                ligands.append(lig_smi)
        except Exception:
            continue
    return ligands


def _analyze_single_ligand(mol):
    """Analyze a single ligand mol for potential donor atoms. Max safety."""
    if mol is None:
        return None
    try:
        atom_counts = Counter()
        donor_est = {}
        for atom in mol.GetAtoms():
            try:
                sym = atom.GetSymbol()
                if sym not in DONOR_ELEMENTS:
                    continue
                atom_counts[sym] += 1
                charge = atom.GetFormalCharge()
                if charge < 0:
                    donor_est[sym] = donor_est.get(sym, 0) + 1
                elif sym in ("N", "P"):
                    try:
                        if atom.GetTotalNumHs() > 0:
                            donor_est[sym] = donor_est.get(sym, 0) + 1
                    except Exception:
                        donor_est[sym] = donor_est.get(sym, 0) + 1
                elif sym in ("Cl", "Br", "I"):
                    if charge < 0:
                        donor_est[sym] = donor_est.get(sym, 0) + 1
            except Exception:
                continue
        if not donor_est:
            donor_est = dict(atom_counts)
        return donor_est if donor_est else None
    except Exception:
        return None


def identify_donor_atoms(smiles, metal, oxidation_state):
    """Identify donor atoms with proper ligand splitting (FIXED version)."""
    try:
        ligand_smiles_list = split_ligands(smiles)
        if not ligand_smiles_list:
            mol = safe_mol_from_smiles(smiles)
            if mol is None:
                return None
            return _analyze_single_ligand(mol)
        merged = Counter()
        for lig_smi in ligand_smiles_list:
            mol = safe_mol_from_smiles(lig_smi)
            if mol is None:
                continue
            donors = _analyze_single_ligand(mol)
            if donors:
                for elem, count in donors.items():
                    merged[elem] += count
        if not merged:
            return None
        donor_est = dict(merged)
        coord_range = COORD_NUM.get((metal, oxidation_state))
        target_coord = coord_range[1] if coord_range else 6
        total = sum(donor_est.values())
        if total > target_coord + 2:
            factor = target_coord / total
            donor_est = {k: max(1, int(v * factor)) for k, v in donor_est.items()}
        return donor_est
    except Exception:
        return None


# ── Worker function for multiprocessing ──

def process_one(row):
    """Process a single complex. Returns (id, donor_json_str) or (id, None)."""
    cid, smiles, metal, ox_state = row
    if not smiles:
        return (cid, None)
    try:
        donors = identify_donor_atoms(smiles, metal, ox_state)
        if donors:
            return (cid, json.dumps(donors))
        return (cid, None)
    except Exception:
        return (cid, None)


# ── Main ──

def main():
    parser = argparse.ArgumentParser(description="Batch donor_atoms enrichment")
    parser.add_argument("--db", default=None, help="Path to biometaldb.sqlite")
    parser.add_argument("--workers", type=int, default=None, help="Number of workers (default: CPU count)")
    parser.add_argument("--batch-size", type=int, default=200, help="Commit batch size")
    parser.add_argument("--dry-run", action="store_true", help="Don't write to DB, just count")
    args = parser.parse_args()

    db_path = args.db or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite"
    )
    n_workers = args.workers or min(cpu_count(), 8)
    batch_size = args.batch_size

    print(f"DB: {db_path}")
    print(f"Workers: {n_workers}")
    print(f"Batch size: {batch_size}")
    print(f"Dry run: {args.dry_run}")
    print()

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Get ALL complexes with SMILES
    cur.execute("""
        SELECT id, smiles_ligands, metal, oxidation_state
        FROM complexes
        WHERE smiles_ligands IS NOT NULL AND smiles_ligands != ''
        ORDER BY id
    """)
    all_rows = cur.fetchall()
    total = len(all_rows)
    print(f"Complexes to process: {total}")

    if args.dry_run:
        print("DRY RUN — no DB changes")
        conn.close()
        return

    # Process in parallel
    t0 = time.time()
    success = 0
    failed = 0
    batch = []

    with Pool(n_workers) as pool:
        for i, (cid, donor_json) in enumerate(pool.imap_unordered(process_one, all_rows, chunksize=50)):
            if donor_json:
                batch.append((donor_json, cid))
                success += 1
            else:
                failed += 1

            # Commit in batches
            if len(batch) >= batch_size:
                cur.executemany("UPDATE complexes SET donor_atoms = ? WHERE id = ?", batch)
                conn.commit()
                batch = []

            # Progress
            if (i + 1) % 500 == 0 or (i + 1) == total:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                eta = (total - i - 1) / rate if rate > 0 else 0
                print(f"  [{i+1}/{total}] ok={success} fail={failed} | {rate:.0f}/s | ETA {eta:.0f}s")
                sys.stdout.flush()

    # Final commit
    if batch:
        cur.executemany("UPDATE complexes SET donor_atoms = ? WHERE id = ?", batch)
        conn.commit()

    elapsed = time.time() - t0

    # Stats
    cur.execute("SELECT COUNT(*) FROM complexes WHERE donor_atoms IS NOT NULL AND donor_atoms != 'None' AND donor_atoms != ''")
    total_filled = cur.fetchone()[0]

    print(f"\n{'='*50}")
    print(f"DONE in {elapsed:.1f}s")
    print(f"Success: {success}")
    print(f"Failed: {failed}")
    print(f"Total with donor_atoms: {total_filled}/{total}")

    # Top patterns
    cur.execute("SELECT donor_atoms, COUNT(*) FROM complexes WHERE donor_atoms IS NOT NULL AND donor_atoms != 'None' GROUP BY donor_atoms ORDER BY COUNT(*) DESC LIMIT 15")
    print(f"\nTop 15 donor_atoms patterns:")
    for da, cnt in cur.fetchall():
        print(f"  {da}: {cnt}")

    conn.close()


if __name__ == "__main__":
    main()
