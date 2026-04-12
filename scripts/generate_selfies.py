#!/usr/bin/env python3
"""Stage 3: Generate SELFIES for all complexes and save mapping."""

import os
import sys
import sqlite3
import csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "patches"))
from selfies_metal import smiles_to_selfies_metal, selfies_to_smiles_metal

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
DB_PATH = os.path.join(DATA_DIR, "biometaldb.sqlite")
SELFIES_CSV = os.path.join(DATA_DIR, "selfies_mapping.csv")


def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Get all complexes with MetalCytoSMILES
    cur.execute("SELECT id, metal_smiles, smiles_ligands, metal, oxidation_state FROM complexes WHERE metal_smiles IS NOT NULL")
    complexes = cur.fetchall()
    print(f"Processing {len(complexes)} complexes...")
    
    success = 0
    fail = 0
    rows = []
    
    for i, (cid, metal_smiles, orig_smiles, metal, ox) in enumerate(complexes):
        if (i + 1) % 1000 == 0:
            print(f"  Progress: {i+1}/{len(complexes)} ({success} ok, {fail} fail)")
        
        try:
            selfies = smiles_to_selfies_metal(metal_smiles)
            roundtrip = selfies_to_smiles_metal(selfies)
            
            # Update DB
            cur.execute("UPDATE complexes SET selfies = ? WHERE id = ?", (selfies, cid))
            
            rows.append({
                'complex_id': cid,
                'metal_smiles': metal_smiles,
                'selfies': selfies,
                'original_smiles': orig_smiles,
                'metal': metal,
                'oxidation_state': ox,
            })
            success += 1
        except Exception as e:
            fail += 1
            if fail <= 5:
                print(f"  FAIL {cid}: {e}")
    
    conn.commit()
    
    # Write CSV mapping
    if rows:
        with open(SELFIES_CSV, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    
    print(f"\n=== Stage 3 Results ===")
    print(f"  Success: {success}")
    print(f"  Failed: {fail}")
    print(f"  SELFIES CSV: {SELFIES_CSV} ({len(rows)} rows)")
    
    # Verify roundtrip for sample
    cur.execute("SELECT metal_smiles, selfies FROM complexes WHERE selfies IS NOT NULL LIMIT 5")
    print(f"\n  Sample roundtrips:")
    for metal_smiles, selfies in cur.fetchall():
        rt = selfies_to_smiles_metal(selfies)
        from rdkit import Chem
        mol1 = Chem.MolFromSmiles(metal_smiles, sanitize=False)
        mol2 = Chem.MolFromSmiles(rt, sanitize=False)
        n1 = mol1.GetNumAtoms() if mol1 else 0
        n2 = mol2.GetNumAtoms() if mol2 else 0
        status = "✓" if n1 == n2 else "✗"
        print(f"    {status} {metal_smiles[:50]}... → {rt[:50]}... (atoms: {n1}→{n2})")
    
    conn.close()
    print(f"\n✓ Stage 3 complete")


if __name__ == "__main__":
    main()
