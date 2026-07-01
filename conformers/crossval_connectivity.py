#!/usr/bin/env python3
"""Connectivity-aware CREST cross-validation. The real question: are CREST's 'conformers'
actually the SAME molecule (same bond graph) as what we serve, or did GFN2 metadynamics
isomerize/fragment the floppy organometallic? Only same-graph CREST frames are valid
comparators. Among valid frames, we then measure ground-state agreement and coverage.

Bond graph = covalent-radius adjacency (0.45 A tolerance), computed on the reference
(our lowest-E opt conformer). A CREST frame is 'same molecule' iff its adjacency matches.
"""
import numpy as np, glob, os, json, sys

SC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/crest_v2_ens"
CD = "/root/conformer-generation/results/conformers"
H2KCAL = 627.5094740631
# covalent radii (A), Cordero-ish, enough for organics + Ir
RCOV = {"H":0.31,"B":0.84,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"P":1.07,"S":1.05,
        "Cl":1.02,"Br":1.20,"I":1.39,"Ir":1.41,"Ru":1.46,"Fe":1.32,"Si":1.11,"Se":1.20}
TOL = 0.45


def frames(p):
    L = open(p).read().splitlines(); F = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); c = L[i+1].strip()
        try: e = float(c.split()[0])
        except Exception: e = None
        A=[];X=[]
        for ln in L[i+2:i+2+n]:
            s=ln.split(); A.append(s[0]); X.append([float(s[1]),float(s[2]),float(s[3])])
        F.append((A,np.array(X),e)); i+=n+2
    return F


def adj(A, X):
    n=len(A); r=np.array([RCOV.get(a,0.77) for a in A])
    D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
    cut=r[:,None]+r[None,:]+TOL
    M=(D<cut)&(D>0.1)
    return M


def kabsch(P,Q):
    P=P-P.mean(0);Q=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(P.T@Q);d=np.sign(np.linalg.det(Wt.T@V.T))
    return float(np.sqrt(((P@(Wt.T@np.diag([1,1,d])@V.T)-Q)**2).sum()/len(P)))


def heavy(A,X):
    m=np.array([a!="H" for a in A]); return X[m]


rows=[]
for f in sorted(glob.glob(os.path.join(SC,"complex_*_crest.xyz"))):
    cid=int(os.path.basename(f).split("_")[1])
    optp=os.path.join(CD,f"complex_{cid}","opt","conformers_opt.xyz")
    if not os.path.exists(optp): continue
    opt=frames(optp); crest=frames(f)
    na=len(opt[0][0])
    refA,refX,_=opt[0]
    refM=adj(refA,refX)                    # bond graph of OUR ground state
    opt_h=[heavy(A,X) for (A,X,e) in opt]
    same=[]; iso=0
    for (A,X,e) in crest:
        if len(A)!=na: continue
        if np.array_equal(adj(A,X),refM): same.append((heavy(A,X),e))
        else: iso+=1
    n_matched=sum(1 for (A,X,e) in crest if len(A)==na)
    emin_same=min((e for (_,e) in same if e is not None), default=None)
    # ground-state agreement: CREST lowest-E same-graph frame vs OUR lowest-E
    gs=None
    if same:
        se=[(e,hX) for (hX,e) in same if e is not None]
        if se:
            _,cr_gs=min(se,key=lambda t:t[0]); gs=round(kabsch(cr_gs,opt_h[0]),2)
    # coverage among SAME-GRAPH low-E (<=3 kcal) CREST frames
    lowE=[(hX,e) for (hX,e) in same if e is not None and (e-emin_same)*H2KCAL<=3.0]
    cov=sum(1 for (hX,e) in lowE if min(kabsch(hX,oX) for oX in opt_h)<=1.0)
    row={"cid":cid,"n_opt":len(opt),"n_crest":len(crest),"n_size_matched":n_matched,
         "same_graph":len(same),"isomerized":iso,
         "iso_pct":round(100*iso/max(1,n_matched),0),
         "groundstate_rmsd_A":gs,"n_samegraph_lowE":len(lowE),"lowE_covered":cov}
    rows.append(row)
    print(f"{cid}: crest_frames={n_matched} | SAME-graph={len(same)} "
          f"isomerized={iso} ({row['iso_pct']:.0f}%) | GS-agree={gs}A | "
          f"low-E same-graph {cov}/{len(lowE)} covered")

json.dump(rows,open(os.path.join(SC,"crossval_connectivity.json"),"w"),indent=1)
if rows:
    tot_iso=sum(r["isomerized"] for r in rows); tot_m=sum(r["n_size_matched"] for r in rows)
    tot_same_low=sum(r["n_samegraph_lowE"] for r in rows); tot_cov=sum(r["lowE_covered"] for r in rows)
    gs_ok=sum(1 for r in rows if r["groundstate_rmsd_A"] is not None and r["groundstate_rmsd_A"]<=1.0)
    print(f"\nSUMMARY: {tot_iso}/{tot_m} CREST frames isomerized ({100*tot_iso/max(1,tot_m):.0f}%) "
          f"-> invalid comparators.")
    print(f"  Ground-state agreement <=1A: {gs_ok}/{len(rows)} complexes.")
    print(f"  Among SAME-GRAPH low-E CREST conformers: {tot_cov}/{tot_same_low} covered by our ensemble.")
