#!/usr/bin/env python3
"""Consolidated CREST cross-validation report (graded thresholds + flexibility).
Verdict inputs per complex:
  - coordination agreement: mean |Ir-donor(ours) - Ir-donor(CREST-min)| over the 6 donors
  - forward coverage: fraction of low-E (<=3 kcal, same bond graph) CREST conformers within
    {1.0, 1.5, 2.0} A heavy-atom of our served ensemble
  - reverse: fraction of OUR conformers within 1.0 A of some CREST conformer
  - flexibility: n_rotatable_bonds (from metadata)
"""
import numpy as np, glob, os, json, sys

SC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/crest_v2_ens"
CD = "/root/conformer-generation/results/conformers"
H2KCAL = 627.5094740631
RCOV = {"H":0.31,"B":0.84,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"P":1.07,"S":1.05,
        "Cl":1.02,"Br":1.20,"I":1.39,"Ir":1.41,"Ru":1.46,"Fe":1.32,"Si":1.11,"Se":1.20}
TOL=0.45

def frames(p):
    L=open(p).read().splitlines();F=[];i=0
    while i<len(L):
        if not L[i].strip(): i+=1; continue
        n=int(L[i].split()[0]);c=L[i+1].strip()
        try:e=float(c.split()[0])
        except Exception:e=None
        A=[];X=[]
        for ln in L[i+2:i+2+n]:
            s=ln.split();A.append(s[0]);X.append([float(x) for x in s[1:4]])
        F.append((A,np.array(X),e));i+=n+2
    return F
def kab(P,Q):
    P=P-P.mean(0);Q=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(P.T@Q);d=np.sign(np.linalg.det(Wt.T@V.T))
    return float(np.sqrt(((P@(Wt.T@np.diag([1,1,d])@V.T)-Q)**2).sum()/len(P)))
def heavy(A,X):
    m=np.array([a!="H" for a in A]);return X[m]
def adj(A,X):
    r=np.array([RCOV.get(a,0.77) for a in A]);D=np.linalg.norm(X[:,None]-X[None],axis=2)
    return (D<(r[:,None]+r[None]+TOL))&(D>0.1)
def ir_donors(A,X):
    ir=A.index("Ir");d=np.linalg.norm(X-X[ir],axis=1);return np.array(sorted(d[(d>0.1)&(d<2.6)]))

rows=[]
for f in sorted(glob.glob(os.path.join(SC,"complex_*_crest.xyz"))):
    cid=int(os.path.basename(f).split("_")[1])
    optp=os.path.join(CD,f"complex_{cid}","opt","conformers_opt.xyz")
    if not os.path.exists(optp):continue
    opt=frames(optp);crest=frames(f);na=len(opt[0][0])
    refA,refX,_=opt[0];refM=adj(refA,refX)
    opt_h=[heavy(A,X) for (A,X,e) in opt]
    same=[(A,X,e) for (A,X,e) in crest if len(A)==na and np.array_equal(adj(A,X),refM)]
    emin=min((e for (_,_,e) in same if e is not None),default=None)
    lowE=[(heavy(A,X),e) for (A,X,e) in same if e is not None and (e-emin)*H2KCAL<=3.0]
    fwd=[min(kab(hX,oX) for oX in opt_h) for (hX,_) in lowE]
    cov={t:(sum(1 for r in fwd if r<=t)/len(fwd) if fwd else None) for t in (1.0,1.5,2.0)}
    rev=[min(kab(oX,heavy(A,X)) for (A,X,e) in same) for oX in opt_h] if same else []
    our_conf=sum(1 for r in rev if r<=1.0)
    # coordination agreement: our GS donors vs CREST-min-E-samegraph donors
    coord_d=None
    se=[(e,A,X) for (A,X,e) in same if e is not None]
    if se:
        _,cA,cX=min(se,key=lambda t:t[0])
        do=ir_donors(refA,refX);dc=ir_donors(cA,cX)
        if len(do)==len(dc) and len(do)>0: coord_d=float(np.abs(do-dc).mean())
    meta={}
    mp=os.path.join(CD,f"complex_{cid}","metadata.json")
    if os.path.exists(mp):
        try:meta=json.load(open(mp))
        except Exception:pass
    row={"cid":cid,"n_rot":meta.get("n_rotatable_bonds"),"n_opt":len(opt),
         "n_crest_samegraph":len(same),"n_lowE":len(lowE),
         "coord_donor_dev_A":round(coord_d,3) if coord_d is not None else None,
         "cov_1.0":round(cov[1.0],2) if cov[1.0] is not None else None,
         "cov_1.5":round(cov[1.5],2) if cov[1.5] is not None else None,
         "cov_2.0":round(cov[2.0],2) if cov[2.0] is not None else None,
         "worst_miss_A":round(max(fwd),2) if fwd else None,
         "our_confirmed":f"{our_conf}/{len(opt)}"}
    rows.append(row)
    print(f"{cid:>5} rot={str(row['n_rot']):>3} opt={row['n_opt']:>3} crest={len(same):>4} "
          f"lowE={len(lowE):>4} | coordΔ={row['coord_donor_dev_A']}A | "
          f"cov@1.0/1.5/2.0={row['cov_1.0']}/{row['cov_1.5']}/{row['cov_2.0']} "
          f"worst={row['worst_miss_A']}A | ours_conf={row['our_confirmed']}")

json.dump(rows,open(os.path.join(SC,"crossval_report.json"),"w"),indent=1)
print("\nwrote", os.path.join(SC,"crossval_report.json"))
