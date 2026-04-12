#!/usr/bin/env python3
"""Stage 2: Generate MetalCytoSMILES, MOL files, and donor atom metadata."""

import os
import sys
import sqlite3
import re
import json
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
DB_PATH = os.path.join(DATA_DIR, "biometaldb.sqlite")
MOL_DIR = os.path.join(DATA_DIR, "mol")
MOL3_DIR = os.path.join(DATA_DIR, "mol3")

# Metal symbols and their atomic numbers
METALS = {
    'Ru': 44, 'Ir': 77, 'Rh': 45, 'Os': 76, 'Re': 75,
    'Fe': 26, 'Pt': 78, 'Pd': 46, 'Cu': 29, 'Zn': 30,
    'Au': 79, 'Ag': 47, 'Co': 27, 'Ni': 28, 'Mn': 25,
    'Cr': 24, 'Mo': 42, 'W': 74, 'V': 23, 'Ti': 22,
}

# Common donor atoms
DONOR_ATOMS = {'N', 'O', 'S', 'P', 'Cl', 'Br', 'I', 'F', 'Se', 'Te', 'As'}


def make_metal_smiles(metal, oxidation_state):
    """Create metal SMILES token like [Ru+2]."""
    if oxidation_state is None:
        return f"[{metal}]"
    ox = int(oxidation_state)
    if ox != 0:
        charge = f"+{ox}" if ox > 0 else str(ox)
        return f"[{metal}{charge}]"
    return f"[{metal}]"


def detect_donor_atoms(smiles_ligands):
    """Detect potential donor atoms in ligand SMILES."""
    donors = defaultdict(int)
    
    # Parse each fragment
    fragments = smiles_ligands.split('.')
    for frag in fragments:
        # Skip counterions
        if re.match(r'^\[(Cl|Br|I|F)-\]$', frag):
            continue
        if re.match(r'^\[.*-\]$', frag) and len(frag) < 10:
            continue
            
        # Find donor atoms
        # Charged atoms: [N+], [O-], [S-], etc.
        for m in re.finditer(r'\[([A-Z][a-z]?)[+-]', frag):
            atom = m.group(1)
            if atom in DONOR_ATOMS:
                donors[atom] += 1
        
        # Neutral N, O, S, P in rings or chains (not in brackets)
        for m in re.finditer(r'(?<!\[)(?<![A-Za-z])([NOSP])(?![A-Za-z\]])(?![+-])', frag):
            donors[m.group(1)] += 1
        
        # Bracketed neutral: [nH], [NH], [OH], etc.
        for m in re.finditer(r'\[([A-Z][a-z]?)[H\d]*\]', frag):
            atom = m.group(1)
            if atom in DONOR_ATOMS:
                donors[atom] += 1
    
    return dict(donors)


def generate_metalcytosmiles(smiles_ligands, metal, oxidation_state, counterion):
    """Generate MetalCytoSMILES: [Metal+charge].ligand1.ligand2.counterion"""
    metal_part = make_metal_smiles(metal, oxidation_state)
    
    # Combine ligands
    parts = [metal_part, smiles_ligands]
    
    # Counterion is already in smiles_ligands in most cases
    # But if counterion is separate, append it
    if counterion and counterion.strip():
        if counterion not in smiles_ligands:
            parts.append(counterion.strip())
    
    return '.'.join(p.strip('.') for p in parts if p.strip())


def generate_mol_file(smiles, output_path, v3000=False):
    """Generate MOL file from MetalCytoSMILES."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "RDKit parse failed"
    
    try:
        AllChem.Compute2DCoords(mol)
    except Exception as e:
        return False, f"2D coords failed: {e}"
    
    try:
        if v3000:
            block = rdmolfiles.MolToV3KMolBlock(mol)
        else:
            block = rdmolfiles.MolToV2KMolBlock(mol)
        
        with open(output_path, 'w') as f:
            f.write(block)
        return True, None
    except Exception as e:
        return False, f"MOL write failed: {e}"


def main():
    os.makedirs(MOL_DIR, exist_ok=True)
    os.makedirs(MOL3_DIR, exist_ok=True)
    
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Get all complexes
    cur.execute("SELECT id, smiles_ligands, metal, oxidation_state, charge_complex FROM complexes")
    complexes = cur.fetchall()
    print(f"Processing {len(complexes)} complexes...")
    
    stats = {'success': 0, 'parse_fail': 0, 'mol_fail': 0, 'mol3_fail': 0}
    donor_stats = defaultdict(int)
    
    for i, (cid, smiles, metal, ox_state, charge) in enumerate(complexes):
        if (i + 1) % 500 == 0:
            print(f"  Progress: {i+1}/{len(complexes)} ({stats['success']} ok, {stats['parse_fail']} fail)")
        
        # Generate MetalCytoSMILES
        metal_smiles = make_metal_smiles(metal, ox_state)
        full_smiles = f"{metal_smiles}.{smiles}"
        
        # Detect donor atoms
        donors = detect_donor_atoms(smiles)
        donors_json = json.dumps(donors) if donors else None
        for atom, count in donors.items():
            donor_stats[atom] += count
        
        # Generate MOL file (V2000)
        mol_path = os.path.join(MOL_DIR, f"complex_{cid}.mol")
        ok, err = generate_mol_file(full_smiles, mol_path, v3000=False)
        
        if not ok and "parse" in err.lower():
            # Try sanitize=False with explicit metal handling
            mol = Chem.MolFromSmiles(full_smiles, sanitize=False)
            if mol is None:
                stats['parse_fail'] += 1
                # Update DB with failure info
                cur.execute("UPDATE complexes SET metal_smiles = ?, donor_atoms = ? WHERE id = ?",
                           (full_smiles, 'PARSE_FAILED', cid))
                continue
        
        if ok:
            # Generate V3000 MOL file
            mol3_path = os.path.join(MOL3_DIR, f"complex_{cid}.mol")
            ok3, _ = generate_mol_file(full_smiles, mol3_path, v3000=True)
            if not ok3:
                stats['mol3_fail'] += 1
            
            stats['success'] += 1
        else:
            stats['mol_fail'] += 1
        
        # Update DB
        cur.execute("UPDATE complexes SET metal_smiles = ?, donor_atoms = ? WHERE id = ?",
                   (full_smiles, donors_json, cid))
    
    conn.commit()
    
    # Report
    print(f"\n=== Stage 2 Results ===")
    print(f"  Success: {stats['success']}")
    print(f"  Parse failed: {stats['parse_fail']}")
    print(f"  MOL write failed: {stats['mol_fail']}")
    print(f"  V3000 failed: {stats['mol3_fail']}")
    print(f"\n  Donor atoms detected:")
    for atom, count in sorted(donor_stats.items(), key=lambda x: -x[1]):
        print(f"    {atom}: {count} occurrences")
    
    # Verify files
    mol_files = [f for f in os.listdir(MOL_DIR) if f.endswith('.mol')]
    mol3_files = [f for f in os.listdir(MOL3_DIR) if f.endswith('.mol')]
    print(f"\n  MOL files: {len(mol_files)}")
    print(f"  MOL3 files: {len(mol3_files)}")
    
    # Spot-check: load a random MOL file back
    if mol_files:
        import random
        sample = random.choice(mol_files)
        sample_path = os.path.join(MOL_DIR, sample)
        mol = Chem.MolFromMolFile(sample_path, sanitize=False)
        if mol:
            print(f"\n  Spot-check {sample}: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() in METALS.values():
                    print(f"    Metal: {atom.GetSymbol()} (charge {atom.GetFormalCharge()})")
        else:
            print(f"\n  Spot-check {sample}: FAILED to reload")
    
    conn.close()
    print(f"\n✓ Stage 2 complete")


if __name__ == "__main__":
    main()
