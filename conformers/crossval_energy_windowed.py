#!/usr/bin/env python3
"""Energy-windowed, heavy-atom cross-validation of our served GFN2-xTB ensemble vs CREST.

Fixes the misleading raw 'cover every CREST conformer' metric:
  - CREST --squick explores a huge HIGH-energy space; we deliberately serve only the
    low-energy window. Demanding we reproduce CREST's high-energy structures is wrong.
  - Use HEAVY-ATOM Kabsch RMSD (drop H noise). Atom indexing is consistent (same input
    topology feeds both Uniconf seeds and CREST), so index-aligned Kabsch is valid.

Per complex we answer the two questions that actually matter for QA:
  Q1 (completeness): does CREST find a LOW-energy conformer (dE<=WIN kcal/mol vs CREST min)
      that is FAR (>1.0 A heavy-atom) from every conformer we serve? -> a real miss.
  Q2 (validity): is each of OUR served conformers confirmed by a nearby CREST conformer?
"""
import numpy as np, glob, os, json, sys

SC  = sys.argv[1] if len(sys.argv) > 1 else "/tmp/crest_v2_ens"
WIN = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0   # low-energy window, kcal/mol
CD  = "/root/conformer-generation/results/conformers"
H2KCAL = 627.5094740631
COV = 1.0   # heavy-atom RMSD (A) below which a conformer is 'covered'


def frames(p):
    L = open(p).read().splitlines(); F = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); comment = L[i + 1].strip()
        try: e = float(comment.split()[0])
        except Exception: e = None
        A = []; X = []
        for ln in L[i + 2:i + 2 + n]:
            s = ln.split(); A.append(s[0]); X.append([float(s[1]), float(s[2]), float(s[3])])
        F.append((A, np.array(X), e)); i += n + 2
    return F


def kabsch(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(P.T @ Q); d = np.sign(np.linalg.det(Wt.T @ V.T))
    return float(np.sqrt(((P @ (Wt.T @ np.diag([1, 1, d]) @ V.T) - Q) ** 2).sum() / len(P)))


def heavy(A, X):
    m = np.array([a != "H" for a in A]); return X[m]


rows = []
for f in sorted(glob.glob(os.path.join(SC, "complex_*_crest.xyz"))):
    cid = int(os.path.basename(f).split("_")[1])
    crest = frames(f)
    optp = os.path.join(CD, f"complex_{cid}", "opt", "conformers_opt.xyz")
    if not os.path.exists(optp):
        print(f"{cid}: no opt ensemble, skip"); continue
    opt = frames(optp)
    na = len(opt[0][0])
    # heavy-atom coords, only frames with matching atom count (same molecule)
    opt_h = [heavy(A, X) for (A, X, e) in opt]
    cr = [(heavy(A, X), e) for (A, X, e) in crest if len(A) == na]
    if not cr:
        print(f"{cid}: no size-matched CREST frames (n_opt_atoms={na})"); continue
    emin = min(e for (_, e) in cr if e is not None)
    # forward: each CREST conformer -> min heavy RMSD to our ensemble, tagged with dE
    fwd = []
    for (cX, e) in cr:
        r = min(kabsch(cX, oX) for oX in opt_h)
        de = (e - emin) * H2KCAL if e is not None else None
        fwd.append((r, de))
    lowE = [(r, de) for (r, de) in fwd if de is not None and de <= WIN]
    n_lowE = len(lowE)
    cov_lowE = sum(1 for (r, de) in lowE if r <= COV)
    miss_lowE = [(round(r, 2), round(de, 2)) for (r, de) in lowE if r > COV]
    worst_lowE_miss = max((r for (r, de) in lowE if r > COV), default=0.0)
    # reverse: each OUR conformer -> nearest CREST conformer (any energy)
    rev = [min(kabsch(oX, cX) for (cX, _) in cr) for oX in opt_h]
    our_confirmed = sum(1 for r in rev if r <= COV)
    row = {
        "cid": cid, "n_opt": len(opt), "n_crest_matched": len(cr),
        "crest_dE_max_kcal": round(max((e - emin) * H2KCAL for (_, e) in cr if e is not None), 1),
        "lowE_window_kcal": WIN, "n_crest_lowE": n_lowE,
        "lowE_covered": cov_lowE, "lowE_missed": n_lowE - cov_lowE,
        "worst_lowE_miss_A": round(worst_lowE_miss, 2),
        "our_confirmed_by_crest": f"{our_confirmed}/{len(opt)}",
        "our_max_rev_rmsd_A": round(max(rev), 2),
    }
    rows.append(row)
    print(f"{cid}: opt={row['n_opt']} crest={row['n_crest_matched']} "
          f"(dE_max={row['crest_dE_max_kcal']}kcal) | lowE(<= {WIN}): "
          f"{cov_lowE}/{n_lowE} covered, worst miss {row['worst_lowE_miss_A']}A | "
          f"OUR confirmed by CREST: {row['our_confirmed_by_crest']} (max {row['our_max_rev_rmsd_A']}A)")

json.dump(rows, open(os.path.join(SC, "crossval_energy_windowed.json"), "w"), indent=1)
# aggregate
if rows:
    fully = sum(1 for r in rows if r["lowE_missed"] == 0)
    tot_low = sum(r["n_crest_lowE"] for r in rows)
    tot_missed = sum(r["lowE_missed"] for r in rows)
    print(f"\nSUMMARY (low-E window <= {WIN} kcal/mol, heavy-atom RMSD, cover<= {COV}A):")
    print(f"  {fully}/{len(rows)} complexes: our served ensemble covers ALL low-E CREST conformers")
    print(f"  {tot_low} low-E CREST conformers total; {tot_missed} not covered by our ensemble")
