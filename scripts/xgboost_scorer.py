#!/usr/bin/env python3
"""
XGBoost coordination mode scorer.

Pipeline:
1. For each complex, use RDMetallics to generate coordination candidates
2. Label candidates using existing donor_atoms as ground truth
3. Extract features per candidate
4. Train XGBoost to distinguish correct from incorrect modes
5. Use model to predict donor atoms for complexes without donor_atoms
"""
import sqlite3
import json
import sys
import os
import numpy as np
import pickle
from collections import Counter

from rdkit import Chem
from rdmetallics.coordinate import coordinate_ligand
from rdmetallics import custom_sanitize, is_transition_metal

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite")
MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "xgboost_scorer.pkl")
TRAINING_DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "xgboost_training.json")

# ==================== Feature tables ====================

ELECTRONEGATIVITY = {
    'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
    'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66,
    'Se': 2.55, 'As': 2.18, 'B': 2.04, 'Si': 1.90,
}

POLARIZABILITY = {
    'H': 0.667, 'C': 1.76, 'N': 1.10, 'O': 0.802, 'F': 0.557,
    'P': 3.63, 'S': 2.90, 'Cl': 2.18, 'Br': 3.05, 'I': 5.35,
    'Se': 3.77, 'As': 4.31,
}

METALS = ['Ru', 'Ir', 'Rh', 'Os', 'Re']

DONOR_ELEMENTS = {'N', 'O', 'S', 'P', 'Cl', 'Br', 'I', 'Se', 'As', 'F'}


# ==================== Data generation ====================

def parse_complex_smiles(smiles):
    """Parse ligand SMILES into individual fragments, skip counterions."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return []
    frags = Chem.GetMolFrags(mol, asMols=True)
    ligands = []
    for frag in frags:
        atoms = [a for a in frag.GetAtoms()]
        symbols = [a.GetSymbol() for a in atoms]
        if len(atoms) == 1 and symbols[0] in ('Cl', 'Br', 'I', 'F'):
            if atoms[0].GetFormalCharge() < 0:
                continue
        if len(atoms) == 1:
            continue
        ligands.append(Chem.MolToSmiles(frag))
    return ligands


def get_coordination_candidates(ligand_smiles, metal_symbol, ox_state, timeout_per_lig=10):
    """Generate coordination candidates using RDMetallics."""
    ligand = Chem.MolFromSmiles(ligand_smiles)
    if ligand is None:
        ligand = Chem.MolFromSmiles(ligand_smiles, sanitize=False)
    if ligand is None:
        return []
    
    ox = ox_state if ox_state is not None else 2
    metal_smi = f"[{metal_symbol}+{ox}]" if ox > 0 else f"[{metal_symbol}{ox}]"
    metal_mol = Chem.MolFromSmiles(metal_smi)
    if metal_mol is None:
        return []
    
    try:
        coord_modes = coordinate_ligand(ligand, metal_mol)
    except:
        return []
    
    # Extract donor sets from each mode
    candidates = []
    seen = set()
    for mode in coord_modes[:30]:
        try:
            custom_sanitize(mode)
        except:
            pass
        
        donors = []
        for atom in mode.GetAtoms():
            if is_transition_metal(atom):
                for n in atom.GetNeighbors():
                    sym = n.GetSymbol()
                    if sym in DONOR_ELEMENTS and not is_transition_metal(n):
                        donors.append(sym)
        
        if donors:
            key = tuple(sorted(donors))
            if key not in seen:
                seen.add(key)
                candidates.append({
                    'donors': sorted(donors),
                    'donor_counts': dict(Counter(donors)),
                    'denticity': len(donors),
                })
    
    return candidates


def extract_features(candidate, metal, ox_state, ligand):
    """
    Extract feature vector for a coordination candidate.
    Returns list of feature values.
    """
    donors = candidate['donors']
    donor_counts = candidate['donor_counts']
    denticity = candidate['denticity']
    
    # Ligand info
    lig_mol = Chem.MolFromSmiles(ligand)
    if lig_mol is None:
        lig_mol = Chem.MolFromSmiles(ligand, sanitize=False)
    
    n_atoms = lig_mol.GetNumAtoms() if lig_mol else 0
    
    # Element counts in ligand
    elem_counts = {}
    if lig_mol:
        for atom in lig_mol.GetAtoms():
            s = atom.GetSymbol()
            elem_counts[s] = elem_counts.get(s, 0) + 1
    
    # Features per donor type (aggregate)
    avg_en = np.mean([ELECTRONEGATIVITY.get(d, 2.0) for d in donors])
    avg_pol = np.mean([POLARIZABILITY.get(d, 2.0) for d in donors])
    
    # Metal one-hot
    metal_features = [1 if m == metal else 0 for m in METALS]
    
    # Donor element one-hot
    donor_elem_features = [donor_counts.get(e, 0) for e in sorted(DONOR_ELEMENTS)]
    
    # Ligand element presence
    ligand_has = [1 if elem_counts.get(e, 0) > 0 else 0 for e in ['N', 'O', 'S', 'P']]
    
    features = (
        metal_features +                          # 5
        [ox_state or 0, denticity, n_atoms] +     # 3
        donor_elem_features +                     # 10
        [avg_en, avg_pol] +                       # 2
        ligand_has +                              # 4
        [len(donors), len(set(donors))]           # 2  (total donors, unique types)
    )
    return features


def label_candidate(candidate, gt_donors):
    """Label: 1 if candidate matches ground truth, 0 otherwise."""
    if not gt_donors:
        return None
    
    cand_counts = candidate['donor_counts']
    gt_counts = gt_donors
    
    # Exact element count match
    if cand_counts == gt_counts:
        return 1
    
    # Partial: same elements but different counts
    cand_elems = set(cand_counts.keys())
    gt_elems = set(gt_counts.keys())
    if cand_elems == gt_elems:
        return 0  # right elements, wrong counts - still negative
    
    return 0


# ==================== Main pipeline ====================

def generate_training_data(n_samples=500, verbose=True):
    """Generate training data from RDMetallics candidates."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    cur.execute("""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, c.donor_atoms
        FROM complexes c
        WHERE c.donor_atoms IS NOT NULL AND c.donor_atoms != ''
        ORDER BY RANDOM()
        LIMIT ?
    """, (n_samples,))
    sample = cur.fetchall()
    conn.close()
    
    if verbose:
        print(f"Processing {len(sample)} complexes...")
    
    X = []
    y = []
    meta = []  # track which complex each row comes from
    stats = {"processed": 0, "with_candidates": 0, "total_candidates": 0, "positive": 0, "negative": 0, "errors": 0}
    
    for i, (cid, smiles, metal, ox_state, donor_atoms_json) in enumerate(sample):
        sys.stdout.write(f"\r[{i+1}/{len(sample)}] candidates:{stats['total_candidates']} pos:{stats['positive']} neg:{stats['negative']}")
        sys.stdout.flush()
        
        # Parse ground truth
        try:
            gt_donors = json.loads(donor_atoms_json) if isinstance(donor_atoms_json, str) else donor_atoms_json
        except:
            stats["errors"] += 1
            continue
        
        # Parse ligands
        ligands = parse_complex_smiles(smiles)
        if not ligands:
            stats["errors"] += 1
            continue
        
        # Generate candidates for each ligand
        all_candidates = []
        for lig in ligands:
            cands = get_coordination_candidates(lig, metal, ox_state)
            for c in cands:
                c['ligand'] = lig
            all_candidates.extend(cands)
        
        if not all_candidates:
            stats["errors"] += 1
            continue
        
        stats["with_candidates"] += 1
        
        # Extract features and label
        for cand in all_candidates:
            label = label_candidate(cand, gt_donors)
            if label is None:
                continue
            
            features = extract_features(cand, metal, ox_state, cand['ligand'])
            X.append(features)
            y.append(label)
            meta.append({"complex_id": cid, "donors": cand['donors'], "denticity": cand['denticity']})
            stats["total_candidates"] += 1
            if label == 1:
                stats["positive"] += 1
            else:
                stats["negative"] += 1
        
        stats["processed"] += 1
    
    if verbose:
        print(f"\n\nDone. Stats: {json.dumps(stats, indent=2)}")
    
    return np.array(X), np.array(y), meta, stats


def train_model(X, y):
    """Train XGBoost model."""
    import xgboost as xgb
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import classification_report, roc_auc_score
    
    # Check class balance
    pos = sum(y == 1)
    neg = sum(y == 0)
    print(f"Class balance: {pos} positive, {neg} negative (ratio: {pos/(pos+neg):.3f})")
    
    if pos < 10:
        print("ERROR: Too few positive samples. Need more training data.")
        return None
    
    # Scale features
    scale_pos_weight = neg / pos if pos > 0 else 1.0
    
    model = xgb.XGBClassifier(
        n_estimators=200,
        max_depth=6,
        learning_rate=0.1,
        scale_pos_weight=scale_pos_weight,
        eval_metric='logloss',
        random_state=42,
        use_label_encoder=False,
    )
    
    # Cross-validation
    if len(X) >= 50:
        scores = cross_val_score(model, X, y, cv=5, scoring='roc_auc')
        print(f"Cross-val AUC: {scores.mean():.3f} ± {scores.std():.3f}")
    
    # Train on all data
    model.fit(X, y)
    
    # Feature importance
    feature_names = (
        [f"metal_{m}" for m in METALS] +
        ["ox_state", "denticity", "n_atoms"] +
        [f"donor_{e}" for e in sorted(DONOR_ELEMENTS)] +
        ["avg_en", "avg_pol"] +
        [f"ligand_has_{e}" for e in ['N', 'O', 'S', 'P']] +
        ["n_donors", "n_unique_types"]
    )
    
    importances = model.feature_importances_
    top_features = sorted(zip(feature_names, importances), key=lambda x: -x[1])[:10]
    print("\nTop 10 features:")
    for name, imp in top_features:
        print(f"  {name}: {imp:.4f}")
    
    return model


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=500, help="Number of complexes to sample")
    parser.add_argument("--skip-training", action="store_true", help="Only generate data, don't train")
    args = parser.parse_args()
    
    # Step 1: Generate training data
    print("=" * 60)
    print("  STEP 1: Generate training data")
    print("=" * 60)
    X, y, meta, stats = generate_training_data(n_samples=args.samples)
    
    if len(X) == 0:
        print("No training data generated. Exiting.")
        return
    
    # Save training data
    training_data = {
        "X": X.tolist(),
        "y": y.tolist(),
        "meta": meta,
        "stats": stats,
    }
    with open(TRAINING_DATA_PATH, "w") as f:
        json.dump(training_data, f)
    print(f"Training data saved: {len(X)} samples")
    
    if args.skip_training:
        print("Skipping training (--skip-training)")
        return
    
    # Step 2: Train model
    print("\n" + "=" * 60)
    print("  STEP 2: Train XGBoost")
    print("=" * 60)
    model = train_model(X, y)
    
    if model is None:
        print("Training failed.")
        return
    
    # Step 3: Save model
    with open(MODEL_PATH, "wb") as f:
        pickle.dump(model, f)
    print(f"\nModel saved to {MODEL_PATH}")
    
    # Step 4: Quick validation - predict for a few complexes
    print("\n" + "=" * 60)
    print("  STEP 3: Quick validation")
    print("=" * 60)
    
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, c.donor_atoms
        FROM complexes c
        WHERE c.donor_atoms IS NOT NULL AND c.donor_atoms != ''
        ORDER BY RANDOM() LIMIT 5
    """)
    
    for cid, smiles, metal, ox_state, gt_json in cur.fetchall():
        gt = json.loads(gt_json)
        ligands = parse_complex_smiles(smiles)
        
        all_cands = []
        for lig in ligands:
            cands = get_coordination_candidates(lig, metal, ox_state)
            for c in cands:
                c['ligand'] = lig
            all_cands.extend(cands)
        
        if not all_cands:
            continue
        
        # Score all candidates
        scored = []
        for cand in all_cands:
            feat = extract_features(cand, metal, ox_state, cand['ligand'])
            prob = model.predict_proba([feat])[0][1]
            scored.append((prob, cand))
        
        scored.sort(key=lambda x: -x[0])
        best = scored[0] if scored else None
        
        print(f"\n#{cid} {metal}({ox_state}) GT={gt}")
        if best:
            print(f"  Predicted: {best[1]['donor_counts']} (score={best[0]:.3f})")
            match = best[1]['donor_counts'] == gt
            print(f"  Match: {'✓' if match else '✗'}")
    
    conn.close()
    print("\nDone!")


if __name__ == "__main__":
    main()
