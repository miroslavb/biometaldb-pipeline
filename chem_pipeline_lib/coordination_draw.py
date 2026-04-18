#!/usr/bin/env python3
"""
Pipeline for drawing coordination compounds.
Supports: MOL block input, SMILES, PubChem lookup.
Outputs: PNG render, MOL V3000, CDXML (ChemDraw).
"""

import os
import json
import subprocess
import requests
from pathlib import Path

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

PIPELINE_DIR = Path(__file__).parent
STRUCTURES_DIR = PIPELINE_DIR / "structures"
RENDERS_DIR = PIPELINE_DIR / "renders"
EXPORTS_DIR = PIPELINE_DIR / "exports"

for d in [STRUCTURES_DIR, RENDERS_DIR, EXPORTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def mol_from_molblock(molblock: str, sanitize: bool = False) -> Chem.Mol:
    """Parse MOL block, skip valence check for metals."""
    mol = Chem.MolFromMolBlock(molblock, sanitize=sanitize, removeHs=False)
    if mol and not sanitize:
        Chem.SanitizeMol(mol, sanitizeOps=(
            Chem.SanitizeFlags.SANITIZE_ALL ^
            Chem.SanitizeFlags.SANITIZE_PROPERTIES
        ))
    return mol


def mol_from_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES with metal-friendly settings."""
    # Try with sanitize=False first for metal-containing molecules
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol:
        try:
            Chem.SanitizeMol(mol, sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_ALL ^
                Chem.SanitizeFlags.SANITIZE_PROPERTIES
            ))
        except Exception:
            pass  # Use unsanitized if sanitization fails
    return mol


def fetch_pubchem_mol(name_or_cid: str) -> dict:
    """Fetch MOL block from PubChem by name or CID."""
    # Check if it's a CID (numeric)
    if name_or_cid.isdigit():
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{name_or_cid}/SDF"
    else:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_or_cid}/SDF"
    
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            molblock = resp.text
            # Also get CID and IUPAC name
            if name_or_cid.isdigit():
                cid = name_or_cid
            else:
                cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_or_cid}/property/IUPACName,IsomericSMILES/JSON"
                cid_resp = requests.get(cid_url, timeout=15)
                if cid_resp.status_code == 200:
                    props = cid_resp.json()['PropertyTable']['Properties'][0]
                    return {
                        'cid': props.get('CID'),
                        'iupac': props.get('IUPACName', ''),
                        'smiles': props.get('IsomericSMILES', ''),
                        'molblock': molblock
                    }
            return {'cid': int(cid), 'molblock': molblock, 'iupac': '', 'smiles': ''}
        else:
            return {'error': f'PubChem returned {resp.status_code}'}
    except Exception as e:
        return {'error': str(e)}


def render_png(mol: Chem.Mol, output_path: str, 
               width: int = 800, height: int = 600,
               legend: str = "", font_size: float = 0.8,
               bond_width: int = 2, highlight_atoms: list = None,
               kekulize: bool = True) -> str:
    """Render molecule to PNG with ACS-style formatting."""
    
    # Generate 2D coordinates
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        rdDepictor.Compute2DCoords(mol)
    
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.bondLineWidth = bond_width
    opts.addAtomIndices = False
    opts.fixedBondLength = 35  # ACS-like spacing
    opts.padding = 0.15
    # Don't show implicit H on C, show NH2 style
    opts.explicitMethyl = False
    opts.includeMetadata = False
    
    if highlight_atoms:
        highlight_bonds = []
        drawer.DrawMolecule(mol, legend=legend,
                          highlightAtoms=highlight_atoms,
                          highlightBonds=highlight_bonds)
    else:
        drawer.DrawMolecule(mol, legend=legend)
    
    drawer.FinishDrawing()
    
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())
    
    return output_path


def mol_to_v3000(mol: Chem.Mol, output_path: str) -> str:
    """Export molecule to MOL V3000 format (full coordination bonds)."""
    molblock = Chem.MolToMolBlock(mol)
    # Convert V2000 to V3000 header
    v3000 = molblock.replace("V2000", "V3000")
    # RDKit's V3000 export is limited, save as V2000 for reliability
    with open(output_path, 'w') as f:
        f.write(molblock)
    return output_path


def mol_to_sdf(mol: Chem.Mol, name: str, output_path: str) -> str:
    """Export to SDF format."""
    writer = Chem.SDWriter(output_path)
    mol.SetProp("_Name", name)
    writer.write(mol)
    writer.close()
    return output_path


def mol_to_cdxml(mol: Chem.Mol, output_path: str) -> str:
    """
    Export to CDXML (ChemDraw) format.
    Generates a minimal CDXML from MOL block.
    """
    molblock = Chem.MolToMolBlock(mol)
    
    # Minimal CDXML template
    cdxml = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">
<CDXML>
<page id="1" BoundingBox="0 0 504 432">
'''
    # Extract atoms and bonds from MOL block
    lines = molblock.strip().split('\n')
    if len(lines) >= 4:
        counts = lines[3].split()
        n_atoms = int(counts[0])
        n_bonds = int(counts[1])
        
        atoms_section = []
        bonds_section = []
        
        for i in range(4, 4 + n_atoms):
            parts = lines[i].split()
            if len(parts) >= 4:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                symbol = parts[3]
                # Scale coordinates for CDXML (roughly 30 units per angstrom)
                cx, cy = 100 + x * 30, 300 - y * 30
                atom_id = i - 3
                atoms_section.append(
                    f'  <n id="{atom_id}" BoundingBox="{cx:.1f} {cy:.1f} {cx+10:.1f} {cy+10:.1f}" '
                    f'p="{cx:.1f} {cy:.1f}"><t><s font="3" size="12" face="96">{symbol}</s></t></n>'
                )
        
        for i in range(4 + n_atoms, 4 + n_atoms + n_bonds):
            parts = lines[i].split()
            if len(parts) >= 3:
                a1, a2 = int(parts[0]), int(parts[1])
                order = int(parts[2])
                bonds_section.append(
                    f'  <b B="{a1}" E="{a2}" Order="{order}"/>'
                )
        
        cdxml += '\n'.join(atoms_section) + '\n'
        cdxml += '\n'.join(bonds_section) + '\n'
    
    cdxml += '''</page>
</CDXML>'''
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(cdxml)
    return output_path


def build_square_planar(metal_symbol: str, ligands: list, 
                         cis_pairs: list = None,
                         metal_charge: int = 0) -> Chem.Mol:
    """
    Build a square planar coordination compound from scratch.
    
    Args:
        metal_symbol: e.g. 'Pt', 'Pd', 'Au'
        ligands: list of (symbol, charge) tuples, e.g. [('Cl', 0), ('Cl', 0), ('N', 0), ('N', 0)]
        cis_pairs: list of tuples indicating which ligands should be cis, e.g. [(0,1)] means ligand 0 and 1 are cis
        metal_charge: formal charge on metal
    
    Returns:
        RDKit Mol object with fixed cis/trans geometry
    """
    from rdkit.Geometry import Point2D
    import math
    
    mol = Chem.RWMol()
    
    # Add metal
    metal_idx = mol.AddAtom(Chem.Atom(metal_symbol))
    mol.GetAtomWithIdx(metal_idx).SetFormalCharge(metal_charge)
    
    # Add ligands
    ligand_indices = []
    for sym, charge in ligands:
        idx = mol.AddAtom(Chem.Atom(sym))
        if charge:
            mol.GetAtomWithIdx(idx).SetFormalCharge(charge)
        ligand_indices.append(idx)
    
    # Add bonds
    for lidx in ligand_indices:
        mol.AddBond(metal_idx, lidx, Chem.BondType.SINGLE)
    
    mol = mol.GetMol()
    
    # Set up coordinates for square planar
    # Default: positions at 0°, 90°, 180°, 270°
    angles = [0, 90, 180, 270]
    
    # If cis_pairs specified, rearrange
    if cis_pairs and len(ligands) == 4:
        # Rearrange so cis pairs are at adjacent positions
        placed = set()
        order = []
        for pair in cis_pairs:
            for p in pair:
                if p not in placed:
                    order.append(p)
                    placed.add(p)
        for i in range(4):
            if i not in placed:
                order.append(i)
        # Reorder ligand indices
        ligand_indices = [ligand_indices[i] for i in order]
    
    coordMap = {metal_idx: Point2D(0.0, 0.0)}
    dist = 1.5
    for i, lidx in enumerate(ligand_indices[:4]):
        angle_rad = math.radians(angles[i])
        x = dist * math.cos(angle_rad)
        y = dist * math.sin(angle_rad)
        coordMap[lidx] = Point2D(x, y)
    
    rdDepictor.Compute2DCoords(mol, coordMap=coordMap)
    return mol


def draw_compound(name: str, 
                  smiles: str = None,
                  molblock: str = None,
                  pubchem_name: str = None,
                  legend: str = None,
                  output_prefix: str = None) -> dict:
    """
    Main entry point: draw a coordination compound.
    
    Provide ONE of: smiles, molblock, pubchem_name.
    Returns dict with paths to generated files.
    """
    if output_prefix is None:
        output_prefix = name.replace(' ', '_').replace('/', '_')
    
    results = {'name': name, 'files': {}}
    
    # Get molecule
    mol = None
    if molblock:
        mol = mol_from_molblock(molblock)
        results['source'] = 'molblock'
    elif smiles:
        mol = mol_from_smiles(smiles)
        results['source'] = 'smiles'
        results['smiles'] = smiles
    elif pubchem_name:
        pc = fetch_pubchem_mol(pubchem_name)
        if 'error' not in pc and 'molblock' in pc:
            mol = mol_from_molblock(pc['molblock'])
            results['source'] = 'pubchem'
            results['pubchem_cid'] = pc.get('cid')
            results['iupac'] = pc.get('iupac', '')
            results['smiles'] = pc.get('smiles', '')
        else:
            results['error'] = pc.get('error', 'PubChem lookup failed')
            return results
    
    if mol is None:
        results['error'] = 'Failed to parse molecule'
        return results
    
    # Save MOL file
    mol_path = str(STRUCTURES_DIR / f"{output_prefix}.mol")
    with open(mol_path, 'w') as f:
        f.write(Chem.MolToMolBlock(mol))
    results['files']['mol'] = mol_path
    
    # Render PNG
    png_path = str(RENDERS_DIR / f"{output_prefix}.png")
    render_png(mol, png_path, legend=legend or name)
    results['files']['png'] = png_path
    
    # Export SDF
    sdf_path = str(EXPORTS_DIR / f"{output_prefix}.sdf")
    mol_to_sdf(mol, name, sdf_path)
    results['files']['sdf'] = sdf_path
    
    # Export CDXML
    cdxml_path = str(EXPORTS_DIR / f"{output_prefix}.cdxml")
    mol_to_cdxml(mol, cdxml_path)
    results['files']['cdxml'] = cdxml_path
    
    results['success'] = True
    results['num_atoms'] = mol.GetNumAtoms()
    results['atoms'] = [(a.GetIdx(), a.GetSymbol(), a.GetFormalCharge()) 
                         for a in mol.GetAtoms()]
    
    return results


# CLI interface
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python coordination_draw.py cisplatin")
        print("  python coordination_draw.py --smiles '[NH3][Pt]([NH3])([Cl])[Cl]' 'Cisplatin'")
        print("  python coordination_draw.py --pubchem 'oxaliplatin'")
        sys.exit(1)
    
    if sys.argv[1] == "--smiles" and len(sys.argv) >= 4:
        smiles = sys.argv[2]
        name = sys.argv[3] if len(sys.argv) > 3 else "compound"
        result = draw_compound(name, smiles=smiles)
    elif sys.argv[1] == "--pubchem" and len(sys.argv) >= 3:
        query = sys.argv[2]
        result = draw_compound(query, pubchem_name=query)
    else:
        query = sys.argv[1]
        result = draw_compound(query, pubchem_name=query)
    
    print(json.dumps(result, indent=2, ensure_ascii=False))
