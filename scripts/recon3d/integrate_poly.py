"""Fold polynuclear-builder structures into the full manifest as recovered.

For each poly_<id>.xyz produced by build_poly, copy it into the viewer's struct/,
attach it as the complex's single 'poly' isomer, flip the record to built, and
clear the old failure verdict. Keeps all existing card metadata (smiles, depiction,
ligand table). Idempotent.

  python integrate_poly.py --out out/full --poly polynuclear/out
"""
from __future__ import annotations
import argparse, glob, json, os, re, shutil
import numpy as np


def _read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    head = L[1]
    syms, P = [], []
    for ln in L[2:2 + n]:
        q = ln.split()
        syms.append(q[0]); P.append([float(q[1]), float(q[2]), float(q[3])])
    chg = mult = None
    m = re.search(r"charge=(-?\d+)", head)
    if m:
        chg = int(m.group(1))
    m = re.search(r"mult=(\d+)", head)
    if m:
        mult = int(m.group(1))
    return n, syms, np.array(P), chg, mult


def _min_heavy(syms, P):
    idx = [i for i, s in enumerate(syms) if s != "H"]
    if len(idx) < 2:
        return 9.9
    A = P[idx]
    d = np.linalg.norm(A[:, None, :] - A[None, :, :], axis=-1)
    np.fill_diagonal(d, 9.9)
    return float(d.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/full")
    ap.add_argument("--poly", default="polynuclear/out")
    args = ap.parse_args()
    man_path = os.path.join(args.out, "manifest.json")
    man = json.load(open(man_path))
    by = {r["id"]: r for r in man["records"]}
    struct_dir = os.path.join(args.out, "struct")
    os.makedirs(struct_dir, exist_ok=True)

    n_int = 0
    for xyz in sorted(glob.glob(os.path.join(args.poly, "poly_*.xyz"))):
        cid = int(re.search(r"poly_(\d+)\.xyz", xyz).group(1))
        if cid not in by:
            continue
        fn = f"poly_{cid}.xyz"
        shutil.copy2(xyz, os.path.join(struct_dir, fn))
        n, syms, P, chg, mult = _read_xyz(xyz)
        clash_ok = _min_heavy(syms, P) > 0.9
        r = by[cid]
        r["status"] = "ok"
        r["isomers"] = [{
            "label": "poly", "valid": bool(clash_ok), "energy": None,
            "total_charge": chg, "mult": mult or 1, "n_atoms": n,
            "method": "GFN2-xTB(constrained)+poly",
            "gates": {"no_clash": bool(clash_ok), "bond_lengths_ok": None},
            "files": {"xyz": fn},
            "note": "polynuclear builder (fragment assembly + harmonic-constrained "
                    "GFN2-xTB); metallocene sandwich held rigid",
        }]
        r["retry_meta"] = {"recovered_by": "polynuclear_builder", "strategies_tried": ["build_poly"]}
        for k in ("fail_reason", "fail_label", "fail_hint", "fail_flags"):
            r.pop(k, None)
        n_int += 1

    # recompute counts
    recs = man["records"]
    import collections
    man["n_ok"] = sum(1 for r in recs if r.get("status") == "ok")
    man["n_isomers_total"] = sum(len(r.get("isomers", [])) for r in recs)
    man["n_valid_struct"] = sum(1 for r in recs for i in r.get("isomers", []) if i.get("valid"))
    man["n_recovered"] = sum(1 for r in recs if (r.get("retry_meta") or {}).get("recovered_by"))
    man["n_still_failed"] = sum(1 for r in recs if r.get("status") != "ok")
    man["n_polynuclear"] = sum(1 for r in recs
                               if (r.get("retry_meta") or {}).get("recovered_by") == "polynuclear_builder")
    man["by_fail_reason"] = dict(collections.Counter(
        r.get("fail_reason") for r in recs if r.get("status") != "ok"))
    json.dump(man, open(man_path, "w"), indent=2)
    print(f"integrated {n_int} polynuclear structures | n_ok={man['n_ok']} "
          f"still_failed={man['n_still_failed']} polynuclear={man['n_polynuclear']}")


if __name__ == "__main__":
    main()
