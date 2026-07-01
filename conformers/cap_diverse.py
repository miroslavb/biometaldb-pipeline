#!/usr/bin/env python3
"""Cap each regenerated conformers.xyz to N diverse conformers via farthest-point
sampling (greedy max-min heavy-atom Kabsch RMSD), always keeping frame 0 (lowest UFF
energy). Preserves the extended/diverse conformers we want rather than just lowest-energy.
Backs up the full ensemble to conformers_full.xyz. Idempotent (re-reads full backup)."""
import sys, os, numpy as np

N = int(os.environ.get("CAP_N", "40"))
CD = "/root/conformer-generation/results/conformers"
IDS = [int(x) for x in sys.argv[1:]] or [int(x) for x in
       "5192 5344 5346 5452 5454 5527 5528 5762 5776 5777 5942 5943 5972 5985 "
       "6058 6110 6111 6112 7153 7678 7812 8525 8944 9122".split()]

def read_frames(p):
    L = open(p).read().splitlines(); F = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); block = L[i:i+n+2]; A = []; X = []
        for ln in block[2:2+n]:
            s = ln.split(); A.append(s[0]); X.append([float(s[1]), float(s[2]), float(s[3])])
        F.append((block[1], A, np.array(X))); i += n + 2
    return F

def kab(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(P.T @ Q); d = np.sign(np.linalg.det(Wt.T @ V.T))
    return float(np.sqrt(((P @ (Wt.T @ np.diag([1,1,d]) @ V.T) - Q)**2).sum()/len(P)))

for cid in IDS:
    d = os.path.join(CD, f"complex_{cid}")
    full = os.path.join(d, "conformers_full.xyz"); cur = os.path.join(d, "conformers.xyz")
    if not os.path.exists(full):
        if not os.path.exists(cur): print(f"{cid}: no conformers.xyz, skip"); continue
        os.rename(cur, full)
    F = read_frames(full)
    if len(F) <= N:
        # nothing to cap; restore full as current
        with open(cur, "w") as f:
            for cm, A, X in F:
                f.write(f"{len(A)}\n{cm}\n")
                for a,(x,y,z) in zip(A,X): f.write(f"{a} {x:.6f} {y:.6f} {z:.6f}\n")
        print(f"{cid}: {len(F)} <= {N}, kept all"); continue
    heavy = [X[np.array([a!='H' for a in A])] for cm,A,X in F]
    sel = [0]                                   # always keep lowest-E
    mind = [kab(heavy[0], heavy[j]) for j in range(len(F))]
    while len(sel) < N:
        nxt = int(np.argmax(mind))
        sel.append(nxt)
        for j in range(len(F)):
            r = kab(heavy[nxt], heavy[j])
            if r < mind[j]: mind[j] = r
        mind[nxt] = -1
    sel = sorted(set(sel))
    with open(cur, "w") as f:
        for k in sel:
            cm, A, X = F[k]
            f.write(f"{len(A)}\n{cm}\n")
            for a,(x,y,z) in zip(A,X): f.write(f"{a} {x:.6f} {y:.6f} {z:.6f}\n")
    print(f"{cid}: {len(F)} -> {len(sel)} diverse (farthest-point)")
