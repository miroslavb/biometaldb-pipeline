#!/usr/bin/env python3
"""CREST + GFN2-xTB cross-validation of the Ir(III) conformer ensembles.
For each complex CREST independently samples conformers; we then measure how well OUR
GFN2-xTB-optimized ensemble covers CREST's conformers (min-RMSD, apples-to-apples), and
also vs the raw UFF ensemble for reference. Resumable per complex. Pure stdlib + numpy.
Configurable via env: CREST_MODE, CREST_SUBSET_FILE, CREST_OUT, CREST_DONE, CREST_SUMMARY.
"""
import os, json, tempfile, subprocess, time, shutil, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

BASE        = os.environ.get("CONFOPT_BASE", "/root/confopt")
CONF_DIR    = os.path.join(BASE, "conformers")
CREST       = os.environ.get("CREST_BIN", "/usr/local/bin/crest")
MODE        = os.environ.get("CREST_MODE", "--quick")
SUBSET      = json.load(open(os.path.join(BASE, os.environ.get("CREST_SUBSET_FILE", "crest_subset.json"))))
CHARGE_MAP  = json.load(open(os.path.join(BASE, "charge_uhf_map.json")))
THREADS     = int(os.environ.get("CREST_THREADS", "12"))
CONCURRENCY = int(os.environ.get("CREST_CONCURRENCY", "4"))
TIMEOUT     = int(os.environ.get("CREST_TIMEOUT", "3600"))
OUT  = os.path.join(BASE, os.environ.get("CREST_OUT", "crest_crossval")); os.makedirs(OUT, exist_ok=True)
DONE = os.path.join(BASE, os.environ.get("CREST_DONE", "CREST_ALL_DONE"))
SUMM = os.path.join(BASE, os.environ.get("CREST_SUMMARY", "crest_crossval_summary.json"))


def read_frames(path):
    L = open(path).read().splitlines(); F = []; i = 0
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


def cover(crest, ref):
    na = len(ref[0][0])
    return [min(kabsch(cX, rX) for (rA, rX) in ref) for (cA, cX) in crest if len(cA) == na]


def summ(c):
    if not c: return {"max": None, "mean": None, "n_uncovered_gt1A": None, "n": 0}
    return {"max": max(c), "mean": float(np.mean(c)), "n_uncovered_gt1A": int(sum(x > 1.0 for x in c)), "n": len(c)}


def run_one(cid):
    outj = os.path.join(OUT, f"complex_{cid}.json")
    if os.path.exists(outj): return {"cid": cid, "status": "skip"}
    try:
        info = CHARGE_MAP[str(cid)]; chrg, uhf = info["chrg"], info["uhf"]
        cdir = os.path.join(CONF_DIR, f"complex_{cid}")
        uni = read_frames(os.path.join(cdir, "conformers.xyz"))
        opt_path = os.path.join(cdir, "opt", "conformers_opt.xyz")
        opt = read_frames(opt_path) if os.path.exists(opt_path) else uni
        sA, sX = uni[0]
        wd = tempfile.mkdtemp(prefix=f"crest_{cid}_", dir="/dev/shm")
        with open(os.path.join(wd, "start.xyz"), "w") as f:
            f.write(f"{len(sA)}\nstart\n")
            for a, (x, y, z) in zip(sA, sX): f.write(f"{a} {x:.6f} {y:.6f} {z:.6f}\n")
        t0 = time.time(); tout = False; rc = None
        try:
            p = subprocess.run([CREST, "start.xyz", MODE, "--chrg", str(chrg), "--uhf", str(uhf),
                                "-T", str(THREADS)], cwd=wd, capture_output=True, text=True, timeout=TIMEOUT)
            rc = p.returncode
            open(os.path.join(OUT, f"complex_{cid}_crest.out"), "w").write(p.stdout[-3000:])
        except subprocess.TimeoutExpired:
            tout = True
        el = time.time() - t0
        cfile = os.path.join(wd, "crest_conformers.xyz")
        crest = read_frames(cfile) if os.path.exists(cfile) else []
        res = {"cid": cid, "status": "ok" if crest else "no_output", "mode": MODE,
               "chrg": chrg, "uhf": uhf, "timeout": tout, "returncode": rc, "elapsed_s": round(el),
               "n_uniconf": len(uni), "n_opt": len(opt), "n_crest": len(crest),
               "cover_vs_opt": summ(cover(crest, opt)), "cover_vs_uff": summ(cover(crest, uni))}
        json.dump(res, open(outj, "w"), indent=1)
        if os.path.exists(cfile): shutil.copy(cfile, os.path.join(OUT, f"complex_{cid}_crest.xyz"))
        shutil.rmtree(wd, ignore_errors=True)
        return res
    except Exception as e:
        return {"cid": cid, "status": "error", "err": repr(e), "tb": traceback.format_exc()[:400]}


def main():
    print(f"[{time.strftime('%H:%M:%S')}] CREST {MODE} on {len(SUBSET)} complexes, "
          f"{CONCURRENCY}x{THREADS}T, timeout={TIMEOUT}s -> {OUT}", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=CONCURRENCY) as ex:
        futs = {ex.submit(run_one, c): c for c in SUBSET}
        for fu in as_completed(futs):
            r = fu.result(); results.append(r); co = r.get("cover_vs_opt", {})
            print(f"  {r['cid']} {r['status']} n_crest={r.get('n_crest')} n_opt={r.get('n_opt')} "
                  f"cover_vs_opt(max={co.get('max')},uncov={co.get('n_uncovered_gt1A')}/{co.get('n')}) "
                  f"{r.get('elapsed_s')}s to={r.get('timeout')}", flush=True)
    json.dump(results, open(SUMM, "w"), indent=1)
    open(DONE, "w").write(json.dumps({"n": len(results)}))
    print("CREST DONE", flush=True)


if __name__ == "__main__":
    main()
