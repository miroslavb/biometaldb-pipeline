#!/usr/bin/env python3
"""Spherical-harmonic topology descriptors for BiometalDB 3D structures.

Two rotation-invariant SH descriptor families are computed per complex:

(A) Coordination-shell Steinhardt order parameters Q_l (l=2,4,6,8) on the
    directions from the metal to its first-shell donor atoms. Q_l is the SH
    power spectrum of the bond-orientation distribution and is the standard
    rotation-invariant fingerprint of local coordination geometry
    (Steinhardt, Nelson & Ronchetti, Phys Rev B 1983).

(B) Whole-molecule SH shape spectrum S_l (l=0..LMAX) of the angular mass
    distribution of all atoms about the molecular centroid -- a global,
    rotation-invariant shape/topology fingerprint in the spirit of the
    SH molecular-surface descriptors (Morris et al. 2005; Mavridis et al. 2007).

Ideal-polyhedron references are computed for validation, each complex is
assigned to its nearest reference by Q-vector distance, and that SH-predicted
geometry is compared with the pipeline's assigned geometry.

Run (in arch-env on T07):
  python sh_topology.py [LIMIT]
Outputs under out/full/sh_analysis/: descriptors.csv, summary.json, 3 PNGs.
"""
import json, os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/root/biometaldb-3d/out/full"
STRUCT = os.path.join(OUT, "struct")
RECORDS = os.path.join(OUT, "records.jsonl")
FIG = os.path.join(OUT, "sh_analysis")
os.makedirs(FIG, exist_ok=True)
LMAX = 10
SHELL_CUT = 2.9          # Angstrom; metal-donor first-shell cutoff
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 0

# ---- real spherical harmonics (scipy >=1.15 sph_harm_y, else legacy sph_harm) ----
try:
    from scipy.special import sph_harm_y
    def Ylm(l, m, polar, azim):
        return sph_harm_y(l, m, polar, azim)          # (n, m, theta=polar, phi=azim)
    _API = "sph_harm_y"
except Exception:
    from scipy.special import sph_harm
    def Ylm(l, m, polar, azim):
        return sph_harm(m, l, azim, polar)            # legacy (m, n, azim, polar)
    _API = "sph_harm"

METALS = {"Ru", "Os", "Ir", "Rh", "Re", "Au", "Pt", "Pd", "Cu", "Fe", "Co", "Ni",
          "Zn", "Mn", "V", "Ti", "Ga", "Gd", "Mo", "Ag", "Cd", "Hg", "Pb", "Sn",
          "Cr", "W", "Ca", "Mg", "Eu", "La", "Tb", "Y", "Zr", "Hf", "Ta", "Sb", "Bi"}


def parse_xyz(path):
    with open(path) as f:
        L = f.read().splitlines()
    n = int(L[0].split()[0])
    syms, xyz = [], []
    for ln in L[2:2 + n]:
        p = ln.split()
        if len(p) >= 4:
            syms.append(p[0]); xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.asarray(xyz, float)


def angles(vec):
    r = np.linalg.norm(vec, axis=1)
    r = np.where(r < 1e-9, 1e-9, r)
    polar = np.arccos(np.clip(vec[:, 2] / r, -1, 1))      # colatitude [0,pi]
    azim = np.arctan2(vec[:, 1], vec[:, 0]) % (2 * np.pi)  # [0,2pi]
    return polar, azim


def Ql(vec, l):
    """Steinhardt order parameter (rotation-invariant SH power) over directions."""
    if len(vec) == 0:
        return float("nan")
    polar, azim = angles(vec)
    qlm = np.array([np.mean(Ylm(l, m, polar, azim)) for m in range(-l, l + 1)])
    return float(np.sqrt(4 * np.pi / (2 * l + 1) * np.sum(np.abs(qlm) ** 2)))


def shape_spectrum(vec, lmax=LMAX):
    polar, azim = angles(vec)
    S = []
    for l in range(lmax + 1):
        qlm = np.array([np.mean(Ylm(l, m, polar, azim)) for m in range(-l, l + 1)])
        S.append(float(4 * np.pi / (2 * l + 1) * np.sum(np.abs(qlm) ** 2)))
    return S


def asphericity(xyz):
    c = xyz - xyz.mean(0)
    g = (c.T @ c) / len(c)
    w = np.sort(np.linalg.eigvalsh(g))[::-1]   # l1>=l2>=l3
    b = w[0] - 0.5 * (w[1] + w[2])             # asphericity
    rg2 = w.sum()
    return float(b / rg2) if rg2 > 1e-9 else 0.0


def ideal_refs():
    s3 = math.sqrt(3)
    R = {
        "linear": np.array([[0, 0, 1.], [0, 0, -1.]]),
        "trigonal": np.array([[1, 0, 0.], [-.5, s3 / 2, 0], [-.5, -s3 / 2, 0]]),
        "tetrahedral": np.array([[1, 1, 1.], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / s3,
        "square_planar": np.array([[1, 0, 0.], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
        "trigonal_bipyramidal": np.array([[0, 0, 1.], [0, 0, -1], [1, 0, 0],
                                          [-.5, s3 / 2, 0], [-.5, -s3 / 2, 0]]),
        "square_pyramidal": np.array([[0, 0, 1.], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
        "octahedral": np.array([[1, 0, 0.], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]),
    }
    return {k: np.array([Ql(v, l) for l in (2, 4, 6, 8)]) for k, v in R.items()}, R


def main():
    refQ, refVec = ideal_refs()
    # ---- self-check: octahedron Q4~0.5092, Q6~0.4406 ----
    oc = refQ["octahedral"]
    print(f"[selfcheck] API={_API} octahedral Q=(Q2,Q4,Q6,Q8)={np.round(oc,4)} "
          f"(canonical simple-cubic/octahedral Q4~0.7637, Q6~0.3536; "
          f"tetrahedral Q4 should be ~0.5092)", flush=True)

    rows, miss = [], 0
    with open(RECORDS) as f:
        recs = [json.loads(l) for l in f]
    ok = [r for r in recs if r.get("status") == "ok" and r.get("isomers")]
    if LIMIT:
        ok = ok[:LIMIT]
    print(f"[run] {len(ok)} ok complexes to analyze", flush=True)

    for r in ok:
        metal = r.get("metal")
        iso = r["isomers"][0]
        xyzf = (iso.get("files") or {}).get("xyz")
        if not xyzf:
            miss += 1; continue
        path = os.path.join(STRUCT, xyzf)
        if not os.path.exists(path):
            miss += 1; continue
        try:
            syms, xyz = parse_xyz(path)
        except Exception:
            miss += 1; continue
        # locate metal atom
        mi = [i for i, s in enumerate(syms) if s == metal]
        if not mi:
            mi = [i for i, s in enumerate(syms) if s in METALS]
        if not mi or len(xyz) < 3:
            miss += 1; continue
        mi = mi[0]
        mpos = xyz[mi]
        d = np.linalg.norm(xyz - mpos, axis=1)
        shell = [i for i in range(len(xyz)) if i != mi and d[i] <= SHELL_CUT and syms[i] != "H"]
        if len(shell) < 2:                       # fallback: nearest non-H atoms (>=CN)
            cn = int(r.get("cn") or 6)
            order = [i for i in np.argsort(d) if i != mi and syms[i] != "H"]
            shell = order[:max(cn, 2)]
        # haptic flag: >=5 carbons close to the metal => eta5-Cp / eta6-arene ring
        n_c_near = sum(1 for i in range(len(xyz)) if i != mi and syms[i] == "C" and d[i] <= 2.5)
        is_haptic = n_c_near >= 5
        cvec = xyz[shell] - mpos
        Qs = {l: Ql(cvec, l) for l in (2, 4, 6, 8)}
        Qvec = np.array([Qs[l] for l in (2, 4, 6, 8)])
        pred = min(refQ, key=lambda k: np.linalg.norm(refQ[k] - Qvec))
        # whole molecule shape spectrum about centroid
        S = shape_spectrum(xyz - xyz.mean(0))
        rows.append({
            "id": r["id"], "metal": metal, "ox": r.get("ox"), "cn": r.get("cn"),
            "geom_assigned": r.get("geometry"), "n_shell": len(shell),
            "is_haptic": int(is_haptic),
            "Q2": Qs[2], "Q4": Qs[4], "Q6": Qs[6], "Q8": Qs[8],
            "sh_pred_geom": pred, "asphericity": asphericity(xyz),
            **{f"S{l}": S[l] for l in range(LMAX + 1)},
        })

    print(f"[run] analyzed {len(rows)} (missing/skipped {miss})", flush=True)
    # ---- CSV ----
    import csv
    cols = list(rows[0].keys())
    with open(os.path.join(FIG, "descriptors.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    metals = [r["metal"] for r in rows]
    import collections
    by_metal = collections.Counter(metals)
    # agreement SH-pred vs assigned (only where assigned maps to a reference name)
    comp = [(r["geom_assigned"], r["sh_pred_geom"]) for r in rows if r["geom_assigned"] in refQ]
    agree = sum(1 for a, b in comp if a == b)
    comp_nh = [(r["geom_assigned"], r["sh_pred_geom"]) for r in rows
               if r["geom_assigned"] in refQ and not r["is_haptic"]]
    agree_nh = sum(1 for a, b in comp_nh if a == b)
    n_haptic = sum(r["is_haptic"] for r in rows)
    # clean coordination sphere: non-haptic AND detected shell size == assigned CN
    clean = [r for r in rows if r["geom_assigned"] in refQ and not r["is_haptic"]
             and r.get("cn") and r["n_shell"] == int(r["cn"])]
    agree_clean = sum(1 for r in clean if r["geom_assigned"] == r["sh_pred_geom"])
    shell_eq_cn = sum(1 for r in rows if r.get("cn") and r["n_shell"] == int(r["cn"]))
    # mean Q by predicted geom
    meanQ = {}
    for g in refQ:
        sub = [r for r in rows if r["sh_pred_geom"] == g]
        if sub:
            meanQ[g] = {q: round(float(np.mean([s[q] for s in sub])), 4)
                        for q in ("Q2", "Q4", "Q6", "Q8")}
    summary = {
        "n_analyzed": len(rows), "by_metal": dict(by_metal.most_common()),
        "ideal_refQ": {k: [round(x, 4) for x in v.tolist()] for k, v in refQ.items()},
        "sh_pred_geom_counts": dict(collections.Counter(r["sh_pred_geom"] for r in rows)),
        "n_haptic": n_haptic, "haptic_frac": round(n_haptic / len(rows), 3),
        "agreement_vs_assigned": {"n_comparable": len(comp), "agree": agree,
                                  "rate": round(agree / len(comp), 3) if comp else None},
        "agreement_nonhaptic": {"n_comparable": len(comp_nh), "agree": agree_nh,
                                "rate": round(agree_nh / len(comp_nh), 3) if comp_nh else None},
        "agreement_clean_sphere": {"n_comparable": len(clean), "agree": agree_clean,
                                   "rate": round(agree_clean / len(clean), 3) if clean else None},
        "mean_n_shell": round(float(np.mean([r["n_shell"] for r in rows])), 2),
        "frac_shell_eq_cn": round(shell_eq_cn / len(rows), 3),
        "mean_Q_by_pred_geom": meanQ,
    }
    json.dump(summary, open(os.path.join(FIG, "summary.json"), "w"), indent=2)
    print("[run] summary:", json.dumps(summary)[:800], flush=True)

    # ===================== FIGURES =====================
    Q4 = np.array([r["Q4"] for r in rows]); Q6 = np.array([r["Q6"] for r in rows])
    # Fig 1: Q4-Q6 map colored by predicted geometry + ideal references
    fig, ax = plt.subplots(figsize=(7.2, 6))
    geos = list(refQ)
    cmap = plt.cm.tab10(np.linspace(0, 1, len(geos)))
    for g, c in zip(geos, cmap):
        idx = [i for i, r in enumerate(rows) if r["sh_pred_geom"] == g]
        if idx:
            ax.scatter(Q4[idx], Q6[idx], s=6, alpha=0.35, color=c, label=f"{g} (n={len(idx)})")
        ax.scatter(refQ[g][1], refQ[g][2], marker="*", s=320, color=c,
                   edgecolor="k", linewidth=1.0, zorder=5)
    ax.set_xlabel("$Q_4$ (coordination shell)"); ax.set_ylabel("$Q_6$ (coordination shell)")
    ax.set_title("BiometalDB coordination topology in SH order-parameter space\n"
                 "(points = complexes, stars = ideal polyhedra)")
    ax.legend(fontsize=7, markerscale=1.5, loc="best")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_Q4Q6.png"), dpi=140); plt.close(fig)

    # Fig 2: mean whole-molecule SH shape spectrum by top metals
    fig, ax = plt.subplots(figsize=(7.2, 5))
    ls = np.arange(1, LMAX + 1)
    for m, _ in by_metal.most_common(6):
        sub = [r for r in rows if r["metal"] == m]
        mean = [np.mean([r[f"S{l}"] for r in sub]) for l in ls]
        ax.plot(ls, mean, marker="o", ms=4, label=f"{m} (n={len(sub)})")
    ax.set_xlabel("degree $l$"); ax.set_ylabel("mean $S_l$ (whole-molecule SH power)")
    ax.set_yscale("log"); ax.set_title("Whole-molecule SH shape spectrum by metal")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_spectrum.png"), dpi=140); plt.close(fig)

    # Fig 3: PCA (numpy SVD) of whole-molecule shape spectra S1..S10, colored by metal
    X = np.array([[r[f"S{l}"] for l in ls] for r in rows])
    X = np.log10(np.clip(X, 1e-8, None))
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    U, sv, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    PC = (X - X.mean(0)) @ Vt[:2].T
    var = (sv[:2] ** 2 / (sv ** 2).sum()) * 100
    fig, ax = plt.subplots(figsize=(7.2, 6))
    for m, _ in by_metal.most_common(6):
        idx = [i for i, r in enumerate(rows) if r["metal"] == m]
        ax.scatter(PC[idx, 0], PC[idx, 1], s=7, alpha=0.4, label=f"{m} (n={len(idx)})")
    ax.set_xlabel(f"PC1 ({var[0]:.0f}%)"); ax.set_ylabel(f"PC2 ({var[1]:.0f}%)")
    ax.set_title("PCA of whole-molecule SH shape spectra (by metal)")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_pca.png"), dpi=140); plt.close(fig)
    print("[run] wrote figures + descriptors.csv + summary.json to", FIG, flush=True)


if __name__ == "__main__":
    main()
