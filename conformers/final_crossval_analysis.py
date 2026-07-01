#!/usr/bin/env python3
"""Final CREST cross-validation analysis with a CONNECTIVITY/COORDINATION filter.
CREST GFN2 metadynamics can isomerize/fragment floppy organometallics, producing
'conformers' that aren't valid conformers of the same molecule (huge RMSD). Before
scoring coverage, keep only CREST conformers whose Ir coordination is intact (6 donors
< 2.5 A) and whose heavy-atom count matches. Coverage = min Kabsch RMSD of each VALID
CREST conformer to our GFN2-xTB-optimized ensemble.
"""
import numpy as np, glob, os, json, sys

SC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/crest_v2_ens"
CD = "/root/conformer-generation/results/conformers"


def frames(p):
    L = open(p).read().splitlines(); F = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); A = []; X = []
        for ln in L[i + 2:i + 2 + n]:
            s = ln.split(); A.append(s[0]); X.append([float(s[1]), float(s[2]), float(s[3])])
        F.append((A, np.array(X))); i += n + 2
    return F


def kabsch(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(P.T @ Q); d = np.sign(np.linalg.det(Wt.T @ V.T))
    return float(np.sqrt(((P @ (Wt.T @ np.diag([1, 1, d]) @ V.T) - Q) ** 2).sum() / len(P)))


def n_ir_donors(A, X):
    if "Ir" not in A: return -1
    ir = A.index("Ir"); d = np.linalg.norm(X - X[ir], axis=1)
    return int(((d > 0.1) & (d < 2.5)).sum())


def summ(c):
    if not c: return {"max": None, "mean": None, "uncov_gt1A": None, "n": 0}
    return {"max": round(max(c), 2), "mean": round(float(np.mean(c)), 2),
            "uncov_gt1A": int(sum(x > 1.0 for x in c)), "n": len(c)}


rows = []
for f in sorted(glob.glob(os.path.join(SC, "complex_*_crest.xyz"))):
    cid = int(os.path.basename(f).split("_")[1])
    crest = frames(f)
    opt = frames(os.path.join(CD, f"complex_{cid}", "opt", "conformers_opt.xyz"))
    na = len(opt[0][0]); ref_don = n_ir_donors(*opt[0])
    valid = [(A, X) for (A, X) in crest if len(A) == na and n_ir_donors(A, X) == ref_don]
    cov_all = [min(kabsch(cX, oX) for (oA, oX) in opt) for (cA, cX) in crest if len(cA) == na]
    cov_val = [min(kabsch(cX, oX) for (oA, oX) in opt) for (cA, cX) in valid]
    rows.append({"cid": cid, "n_opt": len(opt), "n_crest": len(crest),
                 "n_crest_valid": len(valid), "ref_donors": ref_don,
                 "cover_all": summ(cov_all), "cover_valid_only": summ(cov_val)})
    print(f"{cid}: n_opt={len(opt)} n_crest={len(crest)} valid(coord intact)={len(valid)} "
          f"| cover_all={summ(cov_all)} | cover_VALID={summ(cov_val)}")

json.dump(rows, open(os.path.join(SC, "final_crossval.json"), "w"), indent=1)
# aggregate
good = [r for r in rows if r["cover_valid_only"]["n"]]
if good:
    covered = sum(1 for r in good if (r["cover_valid_only"]["uncov_gt1A"] or 0) == 0)
    print(f"\nSUMMARY: {len(rows)} complexes; {sum(r['n_crest_valid'] for r in rows)} valid CREST "
          f"conformers total; {covered}/{len(good)} complexes fully covered (all valid CREST "
          f"conformers <1A from our ensemble).")
