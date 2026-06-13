"""Build a ligand library from the DB: classify each unique ligand by chemical
type (priors + pydentate oracle + SMARTS), count frequency, render 2D, emit JSON.

Usage (arch-env): python ligands_library.py --db pilot.sqlite --oracle /tmp/oracle_full.json --out out/ligands
"""
from __future__ import annotations
import argparse, json, os, sqlite3
from collections import Counter, defaultdict
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
import priors

_NHC = Chem.MolFromSmarts("[#6;X2,X3]([#7])[#7]")          # carbene C between two N
_DIIMINE = Chem.MolFromSmarts("[#7]=[#6]~[#6]=[#7]")        # 1,4-diaza (diimine)
_CARBOXYL = Chem.MolFromSmarts("[CX3](=O)[O-,OX2H1]")
_BETADIKET = Chem.MolFromSmarts("[#8]~[#6]~[#6]~[#6]~[#8]")  # 1,5-O,O (β-diketonate)

GROUP_ORDER = [
    "η6-Arene", "η5-Cyclopentadienyl (Cp/Cp*)", "Carbonyl (CO)",
    "Cyclometalated C^N", "N-heterocyclic carbene (NHC)",
    "N,N-chelate (diimine / polypyridyl)", "Polydentate N-donor (≥3)",
    "N-donor (monodentate)", "N,O-donor (Schiff base / amino-acid)",
    "O,O-chelate (β-diketonate / dicarboxylate)", "O-donor (carboxylate/alkoxide)",
    "Phosphine (monodentate)", "Diphosphine (P,P-chelate)", "P,X-mixed donor",
    "S-donor (thiolate / thioether)", "Halide (coordinating)",
    "Mixed / polydentate", "Solvent / aqua", "Counterion / non-coordinating",
    "Other / unclassified",
]
_SOLVENT = {priors._canon(s) for s in ["O", "CS(C)=O", "CC#N", "CO", "CCO", "C1CCOC1", "CN(C)C=O"]}
_CO = priors._canon("[C-]#[O+]")


def classify(canon, mol, pred):
    if canon in priors.COUNTERIONS:
        return "Counterion / non-coordinating"
    if canon in _SOLVENT:
        return "Solvent / aqua"
    if canon == _CO:
        return "Carbonyl (CO)"
    if mol is not None:
        hap = priors.detect_haptic(mol)
        if hap:
            return "η6-Arene" if hap[0] == "arene" else "η5-Cyclopentadienyl (Cp/Cp*)"
    if canon in priors.HALIDES_CANON:
        return "Halide (coordinating)"
    if not pred or "error" in pred or not pred.get("syms"):
        return "Other / unclassified"
    syms = Counter(pred["syms"])
    nP, nC, nN, nO, nS = syms["P"], syms["C"], syms["N"], syms["O"], syms["S"]
    cn = pred.get("cn", sum(syms.values()))
    if nC >= 1 and mol is not None:
        c_idx = [pred["coord_atoms"][i] for i, s in enumerate(pred["syms"]) if s == "C"]
        is_nhc = any(mol.GetAtomWithIdx(i).GetTotalNumHs() == 0 and
                     sum(1 for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetSymbol() == "N") >= 2
                     for i in c_idx if i < mol.GetNumAtoms())
        if is_nhc or (mol.HasSubstructMatch(_NHC) and nN >= 1):
            return "N-heterocyclic carbene (NHC)"
        if nN >= 1:
            return "Cyclometalated C^N"
    if nP >= 2:
        return "Diphosphine (P,P-chelate)"
    if nP == 1 and (nN + nO + nS) == 0:
        return "Phosphine (monodentate)"
    if nP >= 1:
        return "P,X-mixed donor"
    if nS >= 1 and (nN + nO) == 0:
        return "S-donor (thiolate / thioether)"
    if nN >= 3:
        return "Polydentate N-donor (≥3)"
    if nN == 2 and nO == 0:
        return "N,N-chelate (diimine / polypyridyl)"
    if nN >= 1 and nO >= 1:
        return "N,O-donor (Schiff base / amino-acid)"
    if nO >= 2:
        return "O,O-chelate (β-diketonate / dicarboxylate)"
    if nO == 1 and cn == 1:
        return "O-donor (carboxylate/alkoxide)"
    if nN == 1 and cn == 1:
        return "N-donor (monodentate)"
    if sum(syms.values()) >= 3:
        return "Mixed / polydentate"
    return "Other / unclassified"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--oracle", required=True)
    ap.add_argument("--out", default="out/ligands")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "img"), exist_ok=True)
    oracle = json.load(open(a.oracle))

    def canon(f):
        m = Chem.MolFromSmiles(f, sanitize=True) or Chem.MolFromSmiles(f, sanitize=False)
        return Chem.MolToSmiles(m) if m else f

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    freq, metals, examples = Counter(), defaultdict(set), {}
    for cid, metal, smi in con.execute("SELECT id, metal, smiles_ligands FROM complexes"):
        for frag in set(smi.split(".")):
            frag = frag.strip()
            if not frag:
                continue
            c = canon(frag)
            freq[c] += 1
            metals[c].add(metal)
            examples.setdefault(c, cid)
    con.close()

    ligs = []
    for i, (c, n) in enumerate(sorted(freq.items(), key=lambda kv: -kv[1])):
        mol = Chem.MolFromSmiles(c, sanitize=True) or Chem.MolFromSmiles(c, sanitize=False)
        pred = oracle.get(c)
        grp = classify(c, mol, pred)
        png = f"lig_{i}.png"
        ok = False
        if mol is not None:
            try:
                AllChem.Compute2DCoords(mol)
                Draw.MolToFile(mol, os.path.join(a.out, "img", png), size=(260, 200))
                ok = True
            except Exception:
                ok = False
        ligs.append({
            "idx": i, "smiles": c, "group": grp, "count": n,
            "metals": sorted(metals[c]), "example_id": examples[c],
            "denticity": (pred or {}).get("cn"), "donors": (pred or {}).get("syms"),
            "hemilabile": bool((pred or {}).get("hemilabile")),
            "n_heavy": mol.GetNumHeavyAtoms() if mol else None,
            "png": png if ok else None,
        })
    grp_counts = Counter(l["group"] for l in ligs)
    data = {"n_ligands": len(ligs), "n_complexes": sum(1 for _ in [0]) or 9414,
            "group_order": GROUP_ORDER,
            "group_counts": {g: grp_counts.get(g, 0) for g in GROUP_ORDER if grp_counts.get(g)},
            "ligands": ligs}
    json.dump(data, open(os.path.join(a.out, "ligands.json"), "w"))
    print(f"{len(ligs)} ligands; groups: {dict(data['group_counts'])}")


if __name__ == "__main__":
    main()
