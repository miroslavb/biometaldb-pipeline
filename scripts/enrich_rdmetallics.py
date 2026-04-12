#!/usr/bin/env python3
"""
Enrich BiometalDB donor_atoms using RDMetallics coordination prediction.

Approach:
1. Parse complex SMILES into metal + ligand fragments
2. For each ligand, use RDMetallics to predict coordination modes
3. Extract donor atoms from predicted metal-ligand bonds
4. Compare with existing donor_atoms data
5. Fill gaps where donor_atoms is NULL

D-MPNN validation note: The paper (Moldagulov et al., Angew. Chem. 2026) describes
scoring coordination candidates with D-MPNN. RDMetallics package provides the template
matching step but NOT the D-MPNN scoring. We use frequency-based selection instead:
the most commonly predicted donor atoms across all coordination modes are selected.
"""
import sqlite3
import json
import re
import sys
import os
import signal
from collections import Counter
from functools import wraps
import threading

from rdkit import Chem
from rdmetallics.coordinate import coordinate_ligand
from rdmetallics import custom_sanitize, is_transition_metal

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite")
OUTPUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "rdmetallics_enrichment.json")

# Elements that can be donors
DONOR_ELEMENTS = {'N', 'O', 'S', 'P', 'Se', 'Cl', 'Br', 'I', 'F', 'As', 'Te'}


def parse_complex_smiles(smiles):
    """
    Parse ligand-only SMILES into individual ligand fragments.
    Counterions (Cl-, Br-, etc.) are filtered out.
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
        
        # Skip single-atom counterions (Cl-, Br-, I-, F-)
        if len(atoms) == 1 and symbols[0] in ('Cl', 'Br', 'I', 'F'):
            if atoms[0].GetFormalCharge() < 0:
                continue
        
        # Skip single-atom fragments that aren't meaningful ligands
        if len(atoms) == 1:
            continue
        
        lig_smi = Chem.MolToSmiles(frag)
        if lig_smi:
            ligands.append(lig_smi)
    
    return ligands


def predict_donor_atoms(ligand_smiles, metal_symbol, ox_state=None, max_modes=30):
    """
    Use RDMetallics to predict coordination modes and extract donor atoms.
    Returns dict of {element: count} for predicted donors.
    """
    ligand = Chem.MolFromSmiles(ligand_smiles)
    if ligand is None:
        ligand = Chem.MolFromSmiles(ligand_smiles, sanitize=False)
    if ligand is None:
        return None
    
    # Create metal species molecule
    ox = ox_state if ox_state is not None else 2
    metal_smi = f"[{metal_symbol}+{ox}]" if ox > 0 else f"[{metal_symbol}{ox}]"
    metal_mol = Chem.MolFromSmiles(metal_smi)
    if metal_mol is None:
        return None
    
    try:
        coord_modes = coordinate_ligand(ligand, metal_mol)
    except Exception as e:
        return {"error": str(e)[:200]}
    
    if not coord_modes:
        return None
    
    # Extract donor atoms from all coordination modes
    all_donors = []
    for mode in coord_modes[:max_modes]:
        try:
            custom_sanitize(mode)
        except:
            pass
        
        for atom in mode.GetAtoms():
            if is_transition_metal(atom):
                neighbors = atom.GetNeighbors()
                donors = [n.GetSymbol() for n in neighbors if not is_transition_metal(n)]
                # Filter to donor elements only
                donors = [d for d in donors if d in DONOR_ELEMENTS]
                if donors:
                    all_donors.append(tuple(sorted(donors)))
    
    if not all_donors:
        return None
    
    # Find most common donor set
    donor_counts = Counter(all_donors)
    most_common = donor_counts.most_common(1)[0][0]
    
    # Also get element-level statistics
    element_counts = Counter()
    for donors in all_donors:
        for d in donors:
            element_counts[d] += 1
    
    return {
        "predicted_donors": dict(Counter(most_common)),
        "element_frequency": dict(element_counts.most_common()),
        "n_coordination_modes": len(coord_modes),
        "n_donor_sets": len(donor_counts),
    }


def compare_donor_atoms(existing, predicted):
    """Compare existing vs predicted donor atoms."""
    if not existing or not predicted:
        return "no_comparison"
    
    try:
        existing_dict = json.loads(existing) if isinstance(existing, str) else existing
    except:
        return "parse_error"
    
    predicted_dict = predicted.get("predicted_donors", {})
    
    if not predicted_dict:
        return "no_prediction"
    
    # Normalize: compare element sets
    existing_set = set(existing_dict.keys())
    predicted_set = set(predicted_dict.keys())
    
    if existing_set == predicted_set:
        return "match"
    elif existing_set.issubset(predicted_set):
        return "prediction_has_more"
    elif predicted_set.issubset(existing_set):
        return "existing_has_more"
    elif existing_set & predicted_set:
        return "partial_match"
    else:
        return "mismatch"


def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Get complexes, prioritizing those without donor_atoms
    cur.execute("""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, c.donor_atoms,
               (SELECT COUNT(*) FROM measurements WHERE complex_id = c.id) as n_measurements
        FROM complexes c
        WHERE c.smiles_ligands IS NOT NULL
        ORDER BY 
            CASE WHEN c.donor_atoms IS NULL OR c.donor_atoms = '' THEN 0 ELSE 1 END,
            n_measurements DESC
        LIMIT 100
    """)
    sample = cur.fetchall()
    
    print(f"Processing {len(sample)} complexes with RDMetallics...")
    print("=" * 70)
    
    results = []
    stats = {
        "total": 0,
        "parsed_ok": 0,
        "prediction_ok": 0,
        "match": 0,
        "partial_match": 0,
        "mismatch": 0,
        "prediction_has_more": 0,
        "existing_has_more": 0,
        "no_existing": 0,
        "no_prediction": 0,
        "errors": 0,
    }
    
    for i, (cid, smiles, metal, ox_state, donor_atoms, n_meas) in enumerate(sample):
        stats["total"] += 1
        sys.stdout.write(f"\r[{i+1}/{len(sample)}] Complex #{cid} ({metal}) - OK:{stats['prediction_ok']} ERR:{stats['errors']}")
        sys.stdout.flush()
        
        # Parse ligands from smiles_ligands (no metal in SMILES)
        ligands = parse_complex_smiles(smiles)
        if not ligands:
            stats["errors"] += 1
            continue
        stats["parsed_ok"] += 1
        
        # Use metal and oxidation_state from database
        parsed_metal = metal
        parsed_ox = ox_state
        
        # Predict donors for each ligand
        all_predicted_donors = {}
        ligand_predictions = []
        
        for lig_smi in ligands:
            pred = predict_donor_atoms(lig_smi, parsed_metal, parsed_ox)
            ligand_predictions.append({
                "ligand_smiles": lig_smi[:100],
                "prediction": pred,
            })
            if pred and "predicted_donors" in pred:
                for elem, count in pred["predicted_donors"].items():
                    all_predicted_donors[elem] = all_predicted_donors.get(elem, 0) + count
        
        # Compare with existing
        if donor_atoms:
            comparison = compare_donor_atoms(donor_atoms, {"predicted_donors": all_predicted_donors})
        else:
            comparison = "no_existing"
            stats["no_existing"] += 1
        
        # Update stats
        if comparison == "match":
            stats["match"] += 1
        elif comparison == "partial_match":
            stats["partial_match"] += 1
        elif comparison == "mismatch":
            stats["mismatch"] += 1
        elif comparison == "prediction_has_more":
            stats["prediction_has_more"] += 1
        elif comparison == "existing_has_more":
            stats["existing_has_more"] += 1
        
        if all_predicted_donors:
            stats["prediction_ok"] += 1
        else:
            stats["no_prediction"] += 1
        
        result = {
            "complex_id": cid,
            "metal": parsed_metal,
            "oxidation_state": parsed_ox,
            "n_ligands": len(ligands),
            "existing_donor_atoms": donor_atoms,
            "predicted_donor_atoms": all_predicted_donors if all_predicted_donors else None,
            "comparison": comparison,
            "n_measurements": n_meas,
            "ligand_predictions": ligand_predictions,
        }
        results.append(result)
        
        # Rate limiting (RDKit is local, but still be nice)
        if (i + 1) % 50 == 0:
            conn.commit()
    
    conn.close()
    
    # Save results
    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    
    # Print summary
    print("\n\n" + "=" * 70)
    print("  RDMETALLICS ENRICHMENT SUMMARY")
    print("=" * 70)
    print(f"  Total processed:       {stats['total']}")
    print(f"  Parsed OK:             {stats['parsed_ok']}")
    print(f"  Prediction OK:         {stats['prediction_ok']}")
    print(f"  No prediction:         {stats['no_prediction']}")
    print(f"  Errors:                {stats['errors']}")
    print()
    print(f"  Has existing donors:   {stats['total'] - stats['no_existing']}")
    print(f"  No existing donors:    {stats['no_existing']}")
    print()
    print(f"  === Comparison (where both exist) ===")
    print(f"  Exact match:           {stats['match']}")
    print(f"  Partial match:         {stats['partial_match']}")
    print(f"  Prediction has more:   {stats['prediction_has_more']}")
    print(f"  Existing has more:     {stats['existing_has_more']}")
    print(f"  Mismatch:              {stats['mismatch']}")
    print()
    
    # Show examples of mismatches and fills
    mismatches = [r for r in results if r["comparison"] == "mismatch"]
    no_existing = [r for r in results if r["comparison"] == "no_existing" and r["predicted_donor_atoms"]]
    
    if mismatches[:3]:
        print("=== Examples: mismatches ===")
        for r in mismatches[:3]:
            print(f"  #{r['complex_id']} {r['metal']}({r['oxidation_state']})")
            print(f"    Existing:   {r['existing_donor_atoms']}")
            print(f"    Predicted:  {r['predicted_donor_atoms']}")
    
    if no_existing[:5]:
        print("\n=== Examples: filled gaps (no existing donor_atoms) ===")
        for r in no_existing[:5]:
            print(f"  #{r['complex_id']} {r['metal']}({r['oxidation_state']}): {r['predicted_donor_atoms']}")
    
    print(f"\nResults saved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
