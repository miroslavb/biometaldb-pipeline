"""Integrate simple_builder XYZ structures into the manifest."""
import sys, os, json, glob, re, shutil
import numpy as np
from collections import Counter


def _read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    head = L[1]
    syms, P = [], []
    for ln in L[2:2 + n]:
        q = ln.split()
        syms.append(q[0])
        P.append([float(q[1]), float(q[2]), float(q[3])])
    chg = mult = cn = method = None
    m = re.search(r"charge=(-?\d+)", head)
    if m: chg = int(m.group(1))
    m = re.search(r"mult=(\d+)", head)
    if m: mult = int(m.group(1))
    m = re.search(r"method=(\S+)", head)
    if m: method = m.group(1)
    m = re.search(r"cn=(\d+)", head)
    if m: cn = int(m.group(1))
    return n, syms, np.array(P), chg, mult, method, cn


def _min_heavy(syms, P):
    idx = [i for i, s in enumerate(syms) if s != "H"]
    if len(idx) < 2: return 9.9
    A = P[idx]
    d = np.linalg.norm(A[:, None, :] - A[None, :, :], axis=-1)
    np.fill_diagonal(d, 9.9)
    return float(d.min())


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/root/biometaldb-3d/out/full")
    ap.add_argument("--poly", default="/root/biometaldb-3d/polynuclear/out_simple")
    args = ap.parse_args()

    man_path = os.path.join(args.out, "manifest.json")
    man = json.load(open(man_path))
    by = {r["id"]: r for r in man["records"]}
    struct_dir = os.path.join(args.out, "struct")
    os.makedirs(struct_dir, exist_ok=True)

    n_int = 0
    skipped = 0
    xyz_files = sorted(glob.glob(os.path.join(args.poly, "simple_*.xyz")))
    print(f"Found {len(xyz_files)} simple_*.xyz", flush=True)

    for xyz in xyz_files:
        cid = int(re.search(r"simple_(\d+)\.xyz", xyz).group(1))
        if cid not in by:
            skipped += 1
            continue
        fn = f"simple_{cid}.xyz"
        shutil.copy2(xyz, os.path.join(struct_dir, fn))
        n, syms, P, chg, mult, method, cn = _read_xyz(xyz)
        clash_ok = _min_heavy(syms, P) > 0.9
        r = by[cid]
        r["status"] = "ok"
        r["isomers"] = [{
            "label": "simple",
            "valid": bool(clash_ok),
            "energy": None,
            "total_charge": chg, "mult": mult or 1, "n_atoms": n,
            "method": method or "simple_builder",
            "gates": {"no_clash": bool(clash_ok), "bond_lengths_ok": None},
            "files": {"xyz": fn},
            "note": "simple_builder (RDKit ETKDG + metal-centered assembly, no xTB)",
        }]
        r["retry_meta"] = {"recovered_by": "simple_builder", "strategies_tried": ["build"]}
        for k in ("fail_reason", "fail_label", "fail_hint", "fail_flags"):
            r.pop(k, None)
        n_int += 1

    # Recompute counts
    recs = man["records"]
    man["n_ok"] = sum(1 for r in recs if r.get("status") == "ok")
    man["n_isomers_total"] = sum(len(r.get("isomers", [])) for r in recs)
    man["n_valid_struct"] = sum(1 for r in recs for i in r.get("isomers", []) if i.get("valid"))
    man["n_recovered"] = sum(1 for r in recs if (r.get("retry_meta") or {}).get("recovered_by"))
    man["n_still_failed"] = sum(1 for r in recs if r.get("status") != "ok")
    by_reason = Counter()
    for r in recs:
        if r.get("status") != "ok":
            by_reason[r.get("fail_reason") or "unknown"] += 1
    man["by_fail_reason"] = dict(by_reason)
    json.dump(man, open(man_path, "w"), indent=2)
    print(f"integrated {n_int} simple structures | n_ok={man['n_ok']} "
          f"still_failed={man['n_still_failed']}", flush=True)
    print(f"by_fail_reason: {dict(by_reason)}", flush=True)


if __name__ == "__main__":
    main()
