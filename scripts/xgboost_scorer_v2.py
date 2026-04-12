#!/usr/bin/env python3
"""
XGBoost scorer for RDMetallics coordination candidates.

Approach:
1. For each complex with known donor_atoms:
   a. Parse ligands from smiles_ligands
   b. Run coordinate_ligand on each ligand (with timeout)
   c. Extract donor sets from each coordination mode
   d. Label: 1 if donor set matches donor_atoms, 0 otherwise
2. Extract features for each candidate
3. Train XGBoost to score candidates
4. Predict for complexes without donor_atoms
"""
import sqlite3
import json
import os
import sys
import time
import numpy as np
from collections import Counter
from rdkit import Chem
from rdmetallics.coordinate import coordinate_ligand
from rdmetallics import custom_sanitize, is_transition_metal
import xgboost as xgb
from sklearn.model_selection import cross_val_score
import pickle

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite")
MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "xgboost_scorer_v2.pkl")
DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "xgboost_training_data.json")

# Electronegativity (Pauling)
EN = {'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98, 'P': 2.19,
      'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66, 'Se': 2.55, 'As': 2.18,
      'B': 2.04, 'Si': 1.90, 'Te': 2.10}

# Polarizability (Å³)
POLAR = {'H': 0.67, 'C': 1.76, 'N': 1.10, 'O': 0.80, 'F': 0.56, 'P': 3.63,
         'S': 2.90, 'Cl': 2.18, 'Br': 3.05, 'I': 5.35, 'Se': 3.77, 'As': 4.31,
         'B': 1.60, 'Si': 5.38, 'Te': 5.50}

DONOR_ELEMENTS = {'N', 'O', 'S', 'P', 'Se', 'Cl', 'Br', 'I', 'F', 'As', 'Te'}
METALS = ['Ru', 'Ir', 'Rh', 'Os', 'Re']


def parse_ligands(smiles):
    """Parse ligand SMILES, filter counterions."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return []
    frags = Chem.GetMolFrags(mol, asMols=True)
    ligands = []
    for frag in frags:
        atoms = [a.GetSymbol() for a in frag.GetAtoms()]
        if len(atoms) == 1 and atoms[0] in ('Cl', 'Br', 'I', 'F'):
            continue
        if len(atoms) <= 1:
            continue
        ligands.append(Chem.MolToSmiles(frag))
    return ligands


def donor_set_from_mode(mode_mol):
    """Extract donor element counts from a coordination mode."""
    donors = Counter()
    for atom in mode_mol.GetAtoms():
        if is_transition_metal(atom):
            for nbr in atom.GetNeighbors():
                sym = nbr.GetSymbol()
                if sym in DONOR_ELEMENTS:
                    donors[sym] += 1
    return dict(donors)


def donors_match(predicted, ground_truth):
    """Check if predicted donor set matches ground truth."""
    if not predicted:
        return False
    # Normalize
    p = {k: v for k, v in predicted.items() if v > 0}
    g = {k: v for k, v in ground_truth.items() if v > 0}
    # At least elements match
    return set(p.keys()) == set(g.keys())


def extract_features(ligand_smiles, metal, ox_state, donor_set):
    """Extract features for a coordination candidate."""
    lig = Chem.MolFromSmiles(ligand_smiles)
    if lig is None:
        lig = Chem.MolFromSmiles(ligand_smiles, sanitize=False)
    if lig is None:
        return None

    n_atoms = lig.GetNumAtoms()
    n_bonds = lig.GetNumBonds()
    elements = [a.GetSymbol() for a in lig.GetAtoms()]
    elem_counts = Counter(elements)

    # Ring info
    ring_info = lig.GetRingInfo()
    n_rings = ring_info.NumRings()

    # Ligand flags
    has_N = 1 if 'N' in elem_counts else 0
    has_O = 1 if 'O' in elem_counts else 0
    has_S = 1 if 'S' in elem_counts else 0
    has_P = 1 if 'P' in elem_counts else 0
    has_Cl = 1 if 'Cl' in elem_counts else 0

    # Donor element features
    donor_elem = list(donor_set.keys())[0] if donor_set else 'C'
    d_en = EN.get(donor_elem, 2.5)
    d_polar = POLAR.get(donor_elem, 2.0)

    # Average EN of ligand heavy atoms
    heavy_en = [EN.get(e, 2.5) for e in elements if e != 'H']
    avg_en = np.mean(heavy_en) if heavy_en else 2.5

    # Metal one-hot
    metal_features = {f"metal_{m}": 1 if metal == m else 0 for m in METALS}

    # Denticity
    denticity = sum(donor_set.values())

    # Charge estimate from donor elements
    # O- is common anionic donor, Cl- is anionic
    has_anionic = 1 if any(d in ('O', 'S', 'Cl', 'Br', 'I') for d in donor_set) else 0

    # CN triple bonds (potential nitrile donor)
    has_CN = 1 if lig.HasSubstructMatch(Chem.MolFromSmarts("[#6]#[#7]")) else 0

    # Aromatic N (pyridine-like donor)
    aromatic_N = sum(1 for a in lig.GetAtoms() if a.GetSymbol() == 'N' and a.GetIsAromatic())

    features = {
        'n_atoms': n_atoms,
        'n_bonds': n_bonds,
        'n_rings': n_rings,
        'ox_state': ox_state,
        'denticity': denticity,
        'donor_en': d_en,
        'donor_polar': d_polar,
        'avg_en': avg_en,
        'has_anionic': has_anionic,
        'has_CN': has_CN,
        'aromatic_N': aromatic_N,
        'has_N': has_N,
        'has_O': has_O,
        'has_S': has_S,
        'has_P': has_P,
        'has_Cl': has_Cl,
        'n_N': elem_counts.get('N', 0),
        'n_O': elem_counts.get('O', 0),
        'n_S': elem_counts.get('S', 0),
        'n_P': elem_counts.get('P', 0),
        'donor_N': 1 if 'N' in donor_set else 0,
        'donor_O': 1 if 'O' in donor_set else 0,
        'donor_S': 1 if 'S' in donor_set else 0,
        'donor_P': 1 if 'P' in donor_set else 0,
        'donor_Cl': 1 if 'Cl' in donor_set else 0,
        **metal_features,
    }
    return features


def process_complex(cid, smiles, metal, ox_state, donor_atoms_str):
    """Process one complex: generate candidates, extract features, label."""
    try:
        gt = json.loads(donor_atoms_str) if donor_atoms_str else {}
    except:
        return []

    if not gt:
        return []

    ligands = parse_ligands(smiles)
    if not ligands:
        return []

    metal_mol = Chem.MolFromSmiles(f"[{metal}+{ox_state}]" if ox_state and ox_state > 0 else f"[{metal}+2]")
    if metal_mol is None:
        return []

    rows = []
    for lig_smi in ligands:
        lig = Chem.MolFromSmiles(lig_smi)
        if lig is None:
            continue

        # Run coordinate_ligand with implicit timeout
        try:
            modes = coordinate_ligand(lig, metal_mol)
        except:
            continue

        if not modes:
            continue

        # Collect unique donor sets
        seen = set()
        for mode in modes[:30]:
            try:
                custom_sanitize(mode)
            except:
                pass
            ds = donor_set_from_mode(mode)
            ds_key = tuple(sorted(ds.items()))
            if ds_key in seen or not ds:
                continue
            seen.add(ds_key)

            features = extract_features(lig_smi, metal, ox_state or 2, ds)
            if features is None:
                continue

            label = 1 if donors_match(ds, gt) else 0
            features['label'] = label
            features['complex_id'] = cid
            rows.append(features)

    return rows


def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    # Get complexes with donor_atoms and prefer those with measurements
    cur.execute("""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, c.donor_atoms
        FROM complexes c
        WHERE c.donor_atoms IS NOT NULL AND c.donor_atoms != ''
        ORDER BY c.enrichment_score DESC
        LIMIT 300
    """)
    sample = cur.fetchall()
    conn.close()

    print(f"Processing {len(sample)} complexes for training data...")

    all_rows = []
    skipped = 0
    t0 = time.time()

    for i, (cid, smiles, metal, ox, donors) in enumerate(sample):
        elapsed = time.time() - t0
        if elapsed > 300:  # 5 min total budget
            print(f"\nTimeout at {i}/{len(sample)}")
            break

        sys.stdout.write(f"\r[{i+1}/{len(sample)}] rows={len(all_rows)} skipped={skipped}")
        sys.stdout.flush()

        rows = process_complex(cid, smiles, metal, ox, donors)
        if rows:
            all_rows.extend(rows)
        else:
            skipped += 1

    print(f"\n\nGenerated {len(all_rows)} training rows from {len(sample) - skipped} complexes ({skipped} skipped)")

    if len(all_rows) < 50:
        print("Too few rows for training. Exiting.")
        return

    # Save training data
    with open(DATA_PATH, 'w') as f:
        json.dump(all_rows, f)
    print(f"Training data saved to {DATA_PATH}")

    # Prepare for XGBoost
    feature_names = [k for k in all_rows[0].keys() if k not in ('label', 'complex_id')]
    X = np.array([[r[f] for f in feature_names] for r in all_rows])
    y = np.array([r['label'] for r in all_rows])

    print(f"\nClass balance: positive={y.sum()}/{len(y)} ({y.mean():.1%})")

    # Train
    model = xgb.XGBClassifier(
        n_estimators=100,
        max_depth=4,
        learning_rate=0.1,
        scale_pos_weight=max(1, (y == 0).sum() / max(1, (y == 1).sum())),
        eval_metric='auc',
        random_state=42,
    )

    # Cross-validation
    scores = cross_val_score(model, X, y, cv=5, scoring='roc_auc')
    print(f"Cross-val AUC: {scores.mean():.3f} ± {scores.std():.3f}")

    # Fit on all data
    model.fit(X, y)

    # Feature importance
    importances = model.feature_importances_
    top_idx = np.argsort(importances)[::-1][:15]
    print("\nTop 15 features:")
    for idx in top_idx:
        print(f"  {feature_names[idx]}: {importances[idx]:.4f}")

    # Save model + feature names
    with open(MODEL_PATH, 'wb') as f:
        pickle.dump({'model': model, 'feature_names': feature_names}, f)
    print(f"\nModel saved to {MODEL_PATH}")

    # Validate on a few examples
    print("\n=== Validation ===")
    complex_ids = list(set(r['complex_id'] for r in all_rows))
    for cid in complex_ids[:5]:
        cid_rows = [r for r in all_rows if r['complex_id'] == cid]
        if not cid_rows:
            continue
        X_test = np.array([[r[f] for f in feature_names] for r in cid_rows])
        preds = model.predict_proba(X_test)[:, 1]
        best_idx = preds.argmax()
        gt = [r for r in cid_rows if r['label'] == 1]
        gt_donors = None
        if gt:
            gt_donors = {k: v for k, v in gt[0].items() if k.startswith('donor_') and k != 'donor_en' and k != 'donor_polar' and gt[0][k] == 1}

        best = cid_rows[best_idx]
        pred_donors = {k: v for k, v in best.items() if k.startswith('donor_') and k != 'donor_en' and k != 'donor_polar' and best[k] == 1}

        match = "✓" if best['label'] == 1 else "✗"
        print(f"  #{cid}: predicted={pred_donors} (score={preds[best_idx]:.3f}) {match}")


if __name__ == "__main__":
    main()
