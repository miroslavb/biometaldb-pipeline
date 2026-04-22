#!/usr/bin/env python3
"""Regenerate TUCAN for all complexes in BiometalDB with proper metal-donor bonds."""
import sqlite3, re, sys, os, argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, SanitizeFlags

DB = '/root/.hermes-agent2/biometaldb/data/biometaldb.sqlite'

def parse_ligands(smiles_ligands):
    """Split dot-disconnected SMILES, skip single-atom counterions."""
    parts = smiles_ligands.split('.')
    ligands = []
    for p in parts:
        m = Chem.MolFromSmiles(p, sanitize=False)
        if m and m.GetNumAtoms() > 1:
            ligands.append(p)
    return ligands

def embed_ligand(smi, seed=42):
    """Embed ligand in 3D, return (mol_with_H, num_heavy_atoms)."""
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return None
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = seed
    ok = rdDistGeom.EmbedMolecule(mol, params)
    if ok == -1:
        return None
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except:
        pass
    return mol

def build_mol_and_tucan(ligand_mols, metal_symbol, donor_count):
    """Build MOL V3000 with metal-donor bonds, generate TUCAN."""
    all_atoms = []  # (elem, x, y, z)
    all_bonds = []  # (idx1, idx2)
    frag_sizes = []
    
    for i, mol in enumerate(ligand_mols):
        if mol is None:
            return None, None
        cur_off = len(all_atoms)
        na = mol.GetNumAtoms()
        conf = mol.GetConformer()
        for j in range(na):
            elem = mol.GetAtomWithIdx(j).GetSymbol()
            p = conf.GetAtomPosition(j)
            all_atoms.append((elem, p.x, p.y, p.z))
        for b in mol.GetBonds():
            all_bonds.append((cur_off + b.GetBeginAtomIdx(), cur_off + b.GetEndAtomIdx()))
        frag_sizes.append(na)
    
    if not all_atoms:
        return None, None
    
    # Identify donor N atoms (degree=2, in each fragment)
    donor_global = []
    for i, mol in enumerate(ligand_mols):
        frag_off = sum(frag_sizes[:i])
        found = 0
        for j in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(j)
            if atom.GetSymbol() == 'N' and atom.GetDegree() == 2 and found < 2:
                donor_global.append(frag_off + j)
                found += 1
    
    # If not enough N donors, try O, S, P
    if len(donor_global) < donor_count:
        for i, mol in enumerate(ligand_mols):
            frag_off = sum(frag_sizes[:i])
            for j in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(j)
                if atom.GetSymbol() in ('O', 'S', 'P') and atom.GetDegree() >= 2:
                    if frag_off + j not in donor_global and len(donor_global) < donor_count:
                        donor_global.append(frag_off + j)
    
    # Truncate to expected donor count
    donor_global = donor_global[:donor_count]
    
    if len(donor_global) < min(donor_count, 2):
        return None, None  # not enough donors found
    
    # Place metal at donor centroid
    donor_pos = [np.array(all_atoms[d][1:4]) for d in donor_global if d < len(all_atoms)]
    if len(donor_pos) == 0:
        return None, None  # no valid donor positions
    ctr = np.mean(donor_pos, axis=0)
    if ctr.ndim == 0:  # scalar fallback
        return None, None
    ru_idx = len(all_atoms)
    all_atoms.append((metal_symbol, float(ctr[0]), float(ctr[1]), float(ctr[2])))
    
    # Add metal-donor bonds
    for d in donor_global:
        all_bonds.append((ru_idx, d))
    
    n = len(all_atoms)
    
    # Build MOL V3000
    ml = ["Complex", "  Generated", "",
        f" {n:>3d}{len(all_bonds):>3d}  0  0  0  0  0  0  0  0999 V3000",
        "M  V30 BEGIN CTAB", f"M  V30 COUNTS {n} {len(all_bonds)} 0 0 0", "M  V30 BEGIN ATOM"]
    for i in range(n):
        e, x, y, z = all_atoms[i]
        cs = " CHG=2" if i == ru_idx else ""  # simplified charge
        ml.append(f"M  V30 {i+1} {e} {x:.4f} {y:.4f} {z:.4f} 0{cs}")
    ml.append("M  V30 END ATOM")
    ml.append("M  V30 BEGIN BOND")
    for i, (a, b) in enumerate(all_bonds):
        ml.append(f"M  V30 {i+1} 1 {a+1} {b+1}")
    ml.extend(["M  V30 END BOND", "M  V30 END CTAB", "M  END"])
    mol_text = "\n".join(ml) + "\n"
    
    # Generate TUCAN
    from tucan.io.molfile_reader import graph_from_molfile_text
    from tucan.canonicalization import canonicalize_molecule
    from tucan.serialization import serialize_molecule
    
    try:
        G = graph_from_molfile_text(mol_text)
        Gc = canonicalize_molecule(G)
        tucan_str = serialize_molecule(Gc)
        
        # Verify metal degree
        if G.degree(ru_idx) < 2:
            return None, None
        
        return tucan_str, G.degree(ru_idx)
    except Exception as e:
        return None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', type=int, default=0, help='Skip first N complexes')
    args = parser.parse_args()
    
    conn = sqlite3.connect(DB)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    
    # Get all complexes
    cur.execute("SELECT id, smiles_ligands, metal, oxidation_state, donor_atoms FROM complexes ORDER BY id")
    complexes = cur.fetchall()
    total = len(complexes)
    print(f"Total complexes: {total}, starting from {args.start}")
    
    updated = 0
    failed = 0
    skipped = 0
    
    for idx, cx in enumerate(complexes):
        if idx < args.start:
            skipped += 1
            continue
        if idx % 100 == 0:
            print(f"\n[{idx}/{total}] updated={updated} failed={failed} skipped={skipped}")
        
        cid = cx['id']
        smi = cx['smiles_ligands']
        metal = cx['metal']
        
        if not smi or not metal:
            skipped += 1
            continue
        
        # Parse donor count
        try:
            donors = eval(cx['donor_atoms']) if cx['donor_atoms'] else {}
            donor_count = sum(donors.values())
        except:
            donor_count = 4  # default
        
        # Parse and embed ligands
        lig_smis = parse_ligands(smi)
        if not lig_smis:
            failed += 1
            continue
        
        lig_mols = []
        ok = True
        for i, ls in enumerate(lig_smis):
            mol = embed_ligand(ls, seed=42 + i * 100)
            if mol is None:
                ok = False
                break
            lig_mols.append(mol)
        
        if not ok:
            failed += 1
            continue
        
        # Build and generate TUCAN
        result = build_mol_and_tucan(lig_mols, metal, donor_count)
        if result[0] is None:
            failed += 1
            continue
        
        tucan_str, metal_degree = result
        
        # Update DB
        cur.execute("UPDATE complexes SET tucan = ? WHERE id = ?", (tucan_str, cid))
        conn.commit()
        updated += 1
    
    print(f"\n=== DONE ===")
    print(f"Updated: {updated}/{total}")
    print(f"Failed: {failed}/{total}")
    print(f"Skipped: {skipped}/{total}")
    conn.close()


if __name__ == "__main__":
    main()
