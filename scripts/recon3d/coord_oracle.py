"""Coordination oracle — pydentate ensemble GNN (Kulik group, JACS 2026).

Wraps pydentate's coordination-number / coordinating-atom / hemilability models
to answer "how does this ligand coordinate?" from a SMILES. Used to resolve
denticity/donor-atom ambiguity that SMARTS heuristics get wrong (e.g. curcumin
β-diketonate O^O vs dioxophenanthroline fused ether-O), and to surface hemilabile
ligands + alternative coordination modes (which become review variants).

Models are loaded ONCE (module-level cache) and reused across SMILES.
Output per SMILES: {cn, coord_atoms (indices into the SMILES atom order),
syms, prob, hemilabile (0/1), alternatives}.
"""
from __future__ import annotations
import json
import numpy as np
import torch
from rdkit import Chem

from pydentate import predict as _pred
from pydentate.pydentate_lite import MODELS_DIR

_ARCH = {
    "coordination_number": {"atom_fdim": 133, "bond_fdim": 14, "depth": 6,
                            "dropout": 0.3, "ffn_out_size": 6, "hidden_size": 500},
    "coordinating_atoms": {"atom_fdim": 133, "bond_fdim": 14, "depth": 6,
                           "dropout": 0.35, "ffn_out_size": 1, "hidden_size": 600},
    "hemilability": {"atom_fdim": 133, "bond_fdim": 14, "depth": 6,
                     "dropout": 0.3, "ffn_out_size": 1, "hidden_size": 500},
}
_MODELS = {}
import re as _re


def _load(task):
    state = torch.load(MODELS_DIR / f"{task}.pt", map_location="cpu", weights_only=False)
    loaded = state["state_dict"]
    model = _pred.MoleculeModel(task, _ARCH[task])
    msd = model.state_dict()
    fixed = {}
    for k in loaded:
        if _re.match(r"(encoder\.encoder\.)([Wc])", k):
            nk = k.replace("encoder.encoder", "encoder.encoder.0")
        elif _re.match(r"(^ffn)", k):
            nk = k.replace("ffn", "readout")
        else:
            nk = k
        if nk in msd:
            fixed[nk] = loaded[k]
    msd.update(fixed)
    model.load_state_dict(msd)
    model.eval()
    return model


def _ensure():
    if not _MODELS:
        for t in _ARCH:
            _MODELS[t] = _load(t)
        global _PARAMS
        _PARAMS = _pred.Featurization_parameters()


def predict(smiles: str) -> dict | None:
    """Predict coordination for one ligand SMILES. None if unparseable."""
    _ensure()
    mol = _pred.make_mol(smiles)
    if mol is None:
        return None
    preds = {}
    for task in ("coordination_number", "coordinating_atoms"):
        mg = _pred.MolGraph(mol, _PARAMS)
        with torch.no_grad():
            out = _MODELS[task](mg, None)
        out = out[0] if isinstance(out, list) else out
        preds[task] = out.flatten().tolist()
    cn_probs, atom_probs = preds["coordination_number"], preds["coordinating_atoms"]
    cn_unc = float(np.max([1 - p if p >= 0.5 else p for p in cn_probs]))
    atom_unc = float(np.max([1 - p if p >= 0.5 else p for p in atom_probs]))
    feats = torch.tensor([cn_unc, atom_unc], dtype=torch.float).unsqueeze(0)
    mg = _pred.MolGraph(mol, _PARAMS)
    with torch.no_grad():
        hp = _MODELS["hemilability"](mg, feats)
    hp = (hp[0] if isinstance(hp, list) else hp).flatten().tolist()[0]

    cn = int(np.argmax(cn_probs) + 1)
    atoms = [i for i, a in enumerate(np.round(atom_probs)) if a != 0]
    # internal consistency (mirror pydentate_lite)
    if cn != len(atoms):
        if (cn_unc <= atom_unc) or (len(atoms) > 6):
            atoms = [int(i) for i in np.argsort(atom_probs)[-cn:][::-1]]
        else:
            cn = len(atoms)
    atoms = sorted(int(a) for a in atoms)
    syms = [mol.GetAtomWithIdx(a).GetSymbol() for a in atoms]
    prob = float(np.mean([max(cn_probs)] + [atom_probs[a] for a in atoms])) if atoms else float(max(cn_probs))
    return {"cn": cn, "coord_atoms": atoms, "syms": syms,
            "prob": round(prob, 3), "hemilabile": int(round(hp)),
            "hemi_prob": round(float(hp), 3),
            "cn_unc": round(cn_unc, 3), "atom_unc": round(atom_unc, 3)}


def precompute(smiles_list, out_path=None) -> dict:
    """Batch-predict unique SMILES → {smiles: prediction}. Optionally cache JSON."""
    cache = {}
    if out_path:
        try:
            cache = json.load(open(out_path))
        except Exception:
            cache = {}
    todo = [s for s in dict.fromkeys(smiles_list) if s not in cache]
    for i, s in enumerate(todo):
        try:
            cache[s] = predict(s)
        except Exception as e:  # noqa: BLE001
            cache[s] = {"error": f"{type(e).__name__}: {str(e)[:120]}"}
        if out_path and (i % 25 == 0 or i == len(todo) - 1):
            json.dump(cache, open(out_path, "w"))
    if out_path:
        json.dump(cache, open(out_path, "w"))
    return cache


if __name__ == "__main__":
    import sys
    for s in sys.argv[1:]:
        print(s, "->", predict(s))
