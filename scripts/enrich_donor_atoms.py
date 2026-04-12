#!/usr/bin/env python3
"""
Enrich donor_atoms for complexes that don't have them yet.
Method: SMILES analysis via RDKit to identify potential donor atoms.
"""
import sqlite3
import json
import os
import sys
from collections import Counter

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite")

# Common donor atoms in coordination chemistry
DONOR_ELEMENTS = {"N", "O", "S", "P", "Cl", "Br", "I", "F", "Se", "As"}

# Known coordination number ranges by metal and oxidation state
COORD_NUM = {
    ("Ru", 2): (4, 6),
    ("Ru", 3): (6, 6),
    ("Ir", 3): (6, 6),
    ("Rh", 1): (4, 5),
    ("Rh", 3): (6, 6),
    ("Os", 2): (6, 6),
    ("Os", 3): (6, 6),
    ("Os", 4): (6, 6),
    ("Re", 1): (6, 6),
    ("Re", 2): (6, 6),
    ("Re", 3): (6, 6),
    ("Re", 5): (6, 7),
}


def parse_smiles(smiles):
    """Parse SMILES with RDKit, return mol or None."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
        return mol
    except Exception:
        return None


def split_ligands(smiles):
    """
    Split smiles_ligands into individual ligand fragments.
    Skips single-atom counterions (Cl-, Br-, I-, F-).
    Returns list of ligand SMILES.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return []

    frags = Chem.GetMolFrags(mol, asMols=True)
    ligands = []

    for frag in frags:
        atoms = [a for a in frag.GetAtoms()]
        symbols = [a.GetSymbol() for a in atoms]

        # Skip single-atom counterions
        if len(atoms) == 1 and symbols[0] in ("Cl", "Br", "I", "F"):
            if atoms[0].GetFormalCharge() < 0:
                continue

        # Skip any single-atom fragment
        if len(atoms) == 1:
            continue

        lig_smi = Chem.MolToSmiles(frag)
        if lig_smi:
            ligands.append(lig_smi)

    return ligands


def _analyze_single_ligand(mol):
    """
    Analyze a single ligand mol for potential donor atoms.
    Returns dict {element: count} of estimated donors.
    """
    atom_counts = Counter()
    donor_est = {}

    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in DONOR_ELEMENTS:
            continue

        atom_counts[sym] += 1
        charge = atom.GetFormalCharge()

        # Anionic donors (O-, S-) are very likely coordinating
        if charge < 0:
            donor_est[sym] = donor_est.get(sym, 0) + 1
        # Neutral N with lone pair
        elif sym == "N":
            try:
                if atom.GetTotalNumHs() > 0:
                    donor_est[sym] = donor_est.get(sym, 0) + 1
            except Exception:
                donor_est[sym] = donor_est.get(sym, 0) + 1
        # P with lone pair
        elif sym == "P":
            try:
                if atom.GetTotalNumHs() > 0:
                    donor_est[sym] = donor_est.get(sym, 0) + 1
            except Exception:
                donor_est[sym] = donor_est.get(sym, 0) + 1
        # Halides (Cl-, Br-, I-) with negative charge
        elif sym in ("Cl", "Br", "I"):
            if charge < 0:
                donor_est[sym] = donor_est.get(sym, 0) + 1

    if not donor_est:
        # Fallback: all atoms of donor types
        donor_est = dict(atom_counts)

    return donor_est if donor_est else None


def identify_donor_atoms(smiles, metal, oxidation_state):
    """
    Analyze SMILES to identify potential donor atoms.
    Splits into individual ligands first, analyzes each separately,
    then merges donor counts.
    Returns dict like {"N": 2, "O": 1} or None if can't determine.
    """
    ligand_smiles_list = split_ligands(smiles)

    if not ligand_smiles_list:
        # Fallback: try entire SMILES as one molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        return _analyze_single_ligand(mol)

    # Analyze each ligand separately, merge results
    merged = Counter()
    for lig_smi in ligand_smiles_list:
        mol = Chem.MolFromSmiles(lig_smi)
        if mol is None:
            mol = Chem.MolFromSmiles(lig_smi, sanitize=False)
        if mol is None:
            continue
        donors = _analyze_single_ligand(mol)
        if donors:
            for elem, count in donors.items():
                merged[elem] += count

    if not merged:
        return None

    donor_est = dict(merged)

    # Cap at expected coordination number
    coord_range = COORD_NUM.get((metal, oxidation_state))
    if coord_range:
        target_coord = coord_range[1]
    else:
        target_coord = 6

    total = sum(donor_est.values())
    if total > target_coord + 2:
        factor = target_coord / total
        donor_est = {k: max(1, int(v * factor)) for k, v in donor_est.items()}

    return donor_est


def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Get complexes without donor_atoms
    cur.execute("""
        SELECT id, smiles_ligands, metal, oxidation_state 
        FROM complexes 
        WHERE donor_atoms IS NULL OR donor_atoms = 'None' OR donor_atoms = ''
    """)
    complexes = cur.fetchall()
    print(f"Complexes without donor_atoms: {len(complexes)}")
    
    success = 0
    failed = 0
    none_result = 0
    
    for i, (cid, smiles, metal, ox_state) in enumerate(complexes):
        if i % 500 == 0:
            print(f"\n[{i+1}/{len(complexes)}] success={success} failed={failed} none={none_result}")
            sys.stdout.flush()
        
        if not smiles:
            failed += 1
            continue
        
        n_ligands = len(split_ligands(smiles)) if smiles else 0
        donors = identify_donor_atoms(smiles, metal, ox_state)

        if donors:
            donor_json = json.dumps(donors)
            cur.execute("UPDATE complexes SET donor_atoms = ? WHERE id = ?", (donor_json, cid))
            success += 1
            if i < 5:  # Debug first few
                print(f"  #{cid} {metal}({ox_state}) {n_ligands} ligands -> {donors}")
                sys.stdout.flush()
        else:
            none_result += 1
    
    conn.commit()
    
    # Final stats
    cur.execute("SELECT COUNT(*) FROM complexes WHERE donor_atoms IS NOT NULL AND donor_atoms != 'None' AND donor_atoms != ''")
    total_filled = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM complexes")
    total = cur.fetchone()[0]
    
    print(f"\n=== DONE ===")
    print(f"Success: {success}")
    print(f"Failed: {failed}")
    print(f"No donors found: {none_result}")
    print(f"Total with donor_atoms: {total_filled}/{total} ({total_filled/total*100:.1f}%)")
    
    # Distribution
    cur.execute("SELECT donor_atoms, COUNT(*) FROM complexes WHERE donor_atoms IS NOT NULL AND donor_atoms != 'None' GROUP BY donor_atoms ORDER BY COUNT(*) DESC LIMIT 20")
    print("\n=== Top 20 donor_atoms patterns ===")
    for da, cnt in cur.fetchall():
        print(f"  {da}: {cnt}")
    
    conn.close()


if __name__ == "__main__":
    main()
