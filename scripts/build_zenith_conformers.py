#!/usr/bin/env python3
"""Build zenith_conformers.json — cross-referenced conformer data for the /zenith page."""
import json, os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPT_DIR)  # biometaldb root
CONFORMER_DIR = "/root/conformer-generation/results/conformers"
SUMMARY_PATH = "/root/conformer-generation/results/conformer_summary.json"
XTB_PATH = "/root/conformer-generation/results/xtb_energies.json"
IR_CN_PATH = os.path.join(BASE, "data", "ir_cn_families.json")
OUT_PATH = os.path.join(BASE, "data", "zenith_conformers.json")

# Load conformer summary
with open(SUMMARY_PATH) as f:
    summary = json.load(f)

# Load xtb energies
with open(XTB_PATH) as f:
    xtb_data = json.load(f)
xtb_map = {e['compound_id']: e for e in xtb_data}

# Load ir_cn_families for SMILES/metadata cross-ref
ir_cn_map = {}
if os.path.exists(IR_CN_PATH):
    with open(IR_CN_PATH) as f:
        ir_cn = json.load(f)
    for fk, fam in ir_cn.get('families', {}).items():
        for c in fam.get('complexes', []):
            ir_cn_map[c['id']] = c

# Build merged records
records = []
for r in summary['results']:
    cid = r['compound_id']
    rec = {
        'id': cid,
        'conformer_count': r['conformer_count'],
        'n_rotatable_bonds': r.get('n_rotatable_bonds', 0),
        'status': r['status'],
        'rmsd_max': r.get('rmsd_max', 0),
        'method': r.get('method', ''),
    }
    # Add xtb energies
    if cid in xtb_map:
        e = xtb_map[cid]
        rec['energies_hartree'] = e.get('energies_hartree', [])
        rec['energies_kcal'] = e.get('energies_kcal', [])
        if rec['energies_kcal']:
            rec['energy_min_kcal'] = min(rec['energies_kcal'])
            rec['energy_max_kcal'] = max(rec['energies_kcal'])
            rec['energy_range_kcal'] = rec['energy_max_kcal'] - rec['energy_min_kcal']
    # Add compound metadata
    if cid in ir_cn_map:
        cn = ir_cn_map[cid]
        rec['metal'] = cn.get('metal', '')
        rec['oxidation_state'] = cn.get('oxidation_state', 0)
        rec['charge'] = cn.get('charge', 0)
        rec['mw'] = cn.get('mw', 0)
        rec['smiles_ligands'] = cn.get('smiles_ligands', '')
        rec['min_ic50_dark'] = cn.get('min_ic50_dark')
        rec['cn_ligands'] = cn.get('cn_ligands', [])
        rec['nn_ligand'] = cn.get('nn_ligand', '')
        rec['donor_atoms'] = cn.get('donor_atoms', '')

    # Check if SDF exists
    sdf_path = os.path.join(CONFORMER_DIR, f"complex_{cid}", "conformers_sorted.sdf")
    rec['has_sdf'] = os.path.exists(sdf_path)

    records.append(rec)

# Sort by ID
records.sort(key=lambda x: x['id'])

output = {
    'total_compounds': summary['total_success'],
    'total_conformers': summary['total_conformers'],
    'with_multiple': summary['with_at_least_3_conformers'],
    'rigid_single': summary['rigid_single_conformer'],
    'records': records,
}

with open(OUT_PATH, 'w') as f:
    json.dump(output, f, ensure_ascii=False)

print(f"Built {OUT_PATH}: {len(records)} compounds, {summary['total_conformers']} conformers")
