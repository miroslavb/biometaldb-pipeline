#!/usr/bin/env python3
"""Quality conformer refinement: GFN2-xTB geometry OPTIMIZATION of Uniconf conformers
with the CORRECT per-complex molecular charge/spin, then RMSD+energy dedup and energy
re-ranking. Fixes two defects in the original browser data:
  (1) geometries were UFF (Uniconf), never QM-relaxed;
  (2) single-point energies were computed at charge 0 (these are cationic Ir(III)).
Resumable per complex (opt/DONE sentinel). Pure stdlib + numpy. xtb must be in PATH.
"""
import os, sys, json, glob, tempfile, subprocess, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

BASE      = os.environ.get("CONFOPT_BASE", "/root/confopt")
CONF_DIR  = os.path.join(BASE, "conformers")
XTB       = os.environ.get("XTB_BIN", "/usr/local/bin/xtb")
N_WORKERS = int(os.environ.get("N_WORKERS", "52"))
OPT_LEVEL = os.environ.get("OPT_LEVEL", "normal")
TIMEOUT   = int(os.environ.get("PER_CONF_TIMEOUT", "300"))
SCRATCH   = os.environ.get("SCRATCH", "/dev/shm")
H2KCAL    = 627.5094740631
EDIFF_KCAL, RMSD_THRESH = 0.10, 0.5
CHARGE_MAP = json.load(open(os.path.join(BASE, "charge_uhf_map.json")))


def read_frames(path):
    L = open(path).read().splitlines(); frames = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); blk = L[i:i + n + 2]; ats = []; xs = []
        for ln in blk[2:2 + n]:
            s = ln.split(); ats.append(s[0]); xs.append([float(s[1]), float(s[2]), float(s[3])])
        frames.append((ats, np.array(xs))); i += n + 2
    return frames


def kabsch_rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(Wt.T @ V.T)); D = np.diag([1, 1, d])
    Pr = P @ (Wt.T @ D @ V.T)
    return float(np.sqrt(((Pr - Q) ** 2).sum() / len(P)))


def xtb_opt(ats, X, chrg, uhf, wd):
    with open(os.path.join(wd, "in.xyz"), "w") as f:
        f.write(f"{len(ats)}\nin\n")
        for a, (x, y, z) in zip(ats, X): f.write(f"{a} {x:.6f} {y:.6f} {z:.6f}\n")
    env = os.environ.copy(); env["OMP_NUM_THREADS"] = "1"
    try:
        p = subprocess.run([XTB, "in.xyz", "--gfn", "2", "--opt", OPT_LEVEL,
                            "--chrg", str(chrg), "--uhf", str(uhf)],
                           cwd=wd, capture_output=True, text=True, timeout=TIMEOUT, env=env)
    except subprocess.TimeoutExpired:
        return None
    optf = os.path.join(wd, "xtbopt.xyz")
    if p.returncode != 0 or not os.path.exists(optf): return None
    ats2, X2 = read_frames(optf)[-1]
    cl = open(optf).read().splitlines()[1]; toks = cl.replace(":", " ").split(); E = None
    for i, t in enumerate(toks):
        if t == "energy":
            try: E = float(toks[i + 1])
            except Exception: pass
    if E is None:
        for t in cl.split():
            try: E = float(t); break
            except Exception: pass
    return (ats2, X2, E) if E is not None else None


def ir_donors(ats, X):
    if "Ir" not in ats: return []
    ir = ats.index("Ir"); d = np.linalg.norm(X - X[ir], axis=1)
    return sorted(float(x) for x in d if 0.1 < x < 2.5)


def load_template(sdf):
    m = [x for x in open(sdf).read().split("$$$$") if x.strip()][0].splitlines()
    counts = m[3]; na = int(counts[:3]); nb = int(counts[3:6])
    return counts, na, nb, m[4 + na:4 + na + nb]


def make_sdf(kept, counts, na, bonds):
    out = []
    for ci, (ats, X, E) in enumerate(kept):
        L = [f"conformer_{ci}", "  xtbopt2026        3D", "", counts]
        for a in range(na):
            x, y, z = X[a]
            L.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {ats[a]:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
        L += bonds; L.append("M  END")
        out.append("\n".join(L))
    return "\n$$$$\n".join(out) + "\n$$$$\n"


def process(cid):
    info = CHARGE_MAP.get(str(cid))
    cdir = os.path.join(CONF_DIR, f"complex_{cid}")
    optdir = os.path.join(cdir, "opt"); os.makedirs(optdir, exist_ok=True)
    done = os.path.join(optdir, "DONE")
    if os.path.exists(done): return {"cid": cid, "status": "skip"}
    try:
        chrg, uhf = info["chrg"], info["uhf"]
        frames = read_frames(os.path.join(cdir, "conformers.xyz"))
        opt = []
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            for i, (ats, X) in enumerate(frames):
                wd = os.path.join(td, f"c{i}"); os.makedirs(wd)
                r = xtb_opt(ats, X, chrg, uhf, wd)
                if r: opt.append(r)
        if not opt:
            json.dump({"cid": cid, "error": "all opt failed", "n_in": len(frames)},
                      open(os.path.join(optdir, "opt_meta.json"), "w"))
            open(done, "w").write("failed"); return {"cid": cid, "status": "failed", "n_in": len(frames)}
        opt.sort(key=lambda r: r[2])
        kept = []
        for r in opt:
            if not any(abs(r[2] - k[2]) * H2KCAL < EDIFF_KCAL and kabsch_rmsd(r[1], k[1]) < RMSD_THRESH for k in kept):
                kept.append(r)
        energies = [k[2] for k in kept]; emin = energies[0]
        don = ir_donors(kept[0][0], kept[0][1])
        with open(os.path.join(optdir, "conformers_opt.xyz"), "w") as f:
            for ats, X, E in kept:
                f.write(f"{len(ats)}\nenergy: {E:.8f} Eh\n")
                for a, (x, y, z) in zip(ats, X): f.write(f"{a} {x:.6f} {y:.6f} {z:.6f}\n")
        tmpl = os.path.join(cdir, "conformers.sdf")
        if os.path.exists(tmpl):
            counts, na, nb, bonds = load_template(tmpl)
            open(os.path.join(optdir, "conformers_opt_sorted.sdf"), "w").write(make_sdf(kept, counts, na, bonds))
        meta = {"cid": cid, "chrg": chrg, "uhf": uhf, "n_in": len(frames), "n_opt": len(opt),
                "n_kept": len(kept), "energies_hartree": energies,
                "rel_energies_kcal": [(e - emin) * H2KCAL for e in energies],
                "energy_range_kcal": (max(energies) - min(energies)) * H2KCAL,
                "ir_donors_lowE": don, "n_donors": len(don), "coordination_ok": len(don) == 6}
        json.dump(meta, open(os.path.join(optdir, "opt_meta.json"), "w"), indent=1)
        open(done, "w").write("ok")
        return {"cid": cid, "status": "ok", "n_in": len(frames), "n_kept": len(kept), "n_donors": len(don)}
    except Exception as e:
        return {"cid": cid, "status": "error", "err": repr(e), "tb": traceback.format_exc()[:400]}


def main():
    cids = sorted(int(os.path.basename(d).split("_")[1]) for d in glob.glob(os.path.join(CONF_DIR, "complex_*")))
    print(f"[{time.strftime('%H:%M:%S')}] {len(cids)} complexes, {N_WORKERS} workers, opt={OPT_LEVEL}", flush=True)
    t0 = time.time(); done = 0; agg = {}; kept_tot = 0; badcoord = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(process, c): c for c in cids}
        for fu in as_completed(futs):
            r = fu.result(); done += 1; agg[r["status"]] = agg.get(r["status"], 0) + 1
            if r["status"] == "ok":
                kept_tot += r["n_kept"]
                if r.get("n_donors") != 6: badcoord.append(r["cid"])
            if done % 20 == 0 or r["status"] in ("error", "failed"):
                print(f"[{done}/{len(cids)}] {r['cid']} {r['status']} kept_tot={kept_tot} "
                      f"agg={agg} bad_coord={len(badcoord)} {time.time()-t0:.0f}s", flush=True)
                if r["status"] == "error": print("  ERR", r.get("err"), r.get("tb"), flush=True)
    summary = {"total": len(cids), "agg": agg, "total_kept_conformers": kept_tot,
               "bad_coordination_cids": badcoord, "elapsed_s": round(time.time() - t0),
               "opt_level": OPT_LEVEL, "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
    json.dump(summary, open(os.path.join(BASE, "opt_summary.json"), "w"), indent=1)
    open(os.path.join(BASE, "OPT_ALL_DONE"), "w").write(json.dumps(summary))
    print("DONE", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
