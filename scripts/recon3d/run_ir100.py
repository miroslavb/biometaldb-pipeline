"""Reconstruct the first N iridium complexes -> stereoisomers + TREXIO + manifest.

Pipeline per complex: pydentate oracle -> assign -> build_isomers (cis/trans,
fac/mer) -> validate each -> write XYZ/mol2/TREXIO + 2D PNG -> manifest record.
The manifest feeds review_ui.py (manual-verification interface).

Run (on hive or locally, inside arch-env):
  python run_ir100.py --db pilot.sqlite --limit 100 --workers 16 --out out/ir100
"""
from __future__ import annotations
import argparse, json, os, sqlite3, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from concurrent.futures.process import BrokenProcessPool

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

import assign as A
import build as B
import validate as V
import trexio_writer as TX

OUT = "out/ir100"


def _canon(frag):
    m = Chem.MolFromSmiles(frag, sanitize=True) or Chem.MolFromSmiles(frag, sanitize=False)
    return Chem.MolToSmiles(m) if m else frag


def render_png(cid, metal, ox, smiles_ligands, path):
    try:
        smi = f"[{metal}+{ox}]." + smiles_ligands if ox else f"[{metal}]." + smiles_ligands
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            m = Chem.MolFromSmiles(smiles_ligands, sanitize=False)
        AllChem.Compute2DCoords(m)
        Draw.MolToFile(m, path, size=(520, 380))
        return os.path.basename(path)
    except Exception:
        return None


def work(args):
    cid, metal, ox, charge, smiles, donor_json, oracle, outdir = args
    rec = {"id": cid, "metal": metal, "ox": ox, "charge_complex": charge,
           "smiles_ligands": smiles, "isomers": []}
    try:
        donor_atoms = json.loads(donor_json) if donor_json else None
        asg = A.assign(metal, ox, smiles, donor_atoms, oracle=oracle)
        rec.update({"geometry": asg.geometry, "cn": asg.final_cn, "pref_cn": asg.pref_cn,
                    "spin": asg.spin, "assign_conf": asg.confidence, "flags": asg.flags,
                    "hemilabile": asg.hemilabile, "needs_xtb": asg.needs_xtb,
                    "ligands": [{"smiles": u.smiles, "role": u.role, "sites": u.sites,
                                 "coordList": u.pocket, "comp": dict(u.pocket_comp),
                                 "via": ("oracle" if u.oracle else u.role),
                                 "hemilabile": bool(u.oracle and u.oracle.get("hemilabile")),
                                 "oracle": u.oracle} for u in asg.units],
                    "spectators": asg.spectators})
        rec["png"] = render_png(cid, metal, ox, smiles, os.path.join(outdir, "img", f"c{cid}.png"))
        t = time.time()
        isos = B.build_isomers(asg, mode="xtb", max_isomers=2)
        rec["build_s"] = round(time.time() - t, 1)
        for k, res in enumerate(isos):
            label = res.get("isomer", f"iso{k+1}")
            passed, gates, neigh = V.validate(res, asg)
            stem = os.path.join(outdir, "struct", f"c{cid}_{label}")
            paths = B.write_outputs(res, asg, os.path.join(outdir, "struct"), f"{cid}_{label}")
            tx = None
            try:
                tx = TX.write_trexio(stem + ".h5", res["symbols"], res["positions"],
                                     res.get("total_charge"), res.get("n_unpaired"),
                                     title=f"BiometalDB #{cid} {metal}({ox}) {label}")
            except Exception as e:  # noqa: BLE001
                rec.setdefault("trexio_errors", []).append(f"{label}: {type(e).__name__} {str(e)[:80]}")
            rec["isomers"].append({
                "label": label, "valid": passed, "energy": res.get("energy"),
                "total_charge": res.get("total_charge"), "mult": (res.get("n_unpaired") or 0) + 1,
                "n_atoms": res.get("n_atoms"), "method": res.get("method"),
                "stereo_sig": res.get("stereo_sig"),
                "metal_neighbors": [[s, d] for _, s, d in neigh],
                "gates": {k2: gates[k2] for k2 in ("no_clash", "bond_lengths_ok",
                          "min_heavy_dist", "cn_got", "comp_got") if k2 in gates},
                "files": {**{k2: os.path.basename(v) for k2, v in paths.items()},
                          **({"trexio": os.path.basename(tx)} if tx else {})},
            })
        rec["status"] = "ok" if rec["isomers"] else "no_structure"
    except Exception as e:  # noqa: BLE001
        rec["status"] = "exception"
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-600:]
    return rec


def _kill_pool(ex):
    """Force-kill all worker processes of a (possibly hung) pool. shutdown() alone
    will NOT terminate a worker stuck in native code (xtb/Architector), so we kill
    the OS processes directly before tearing the executor down."""
    for proc in list(getattr(ex, "_processes", {}).values()):
        try:
            proc.kill()
        except Exception:  # noqa: BLE001
            pass


def _fast_pass(remaining, workers, timeout, sink, drop):
    """Sliding-window pool over `remaining` (≤ `workers` tasks outstanding, so one
    bad complex can't implicate thousands). Each completed record is sink()'d and
    drop()'d immediately. Returns (suspects, problem):
      problem = None  -> drained cleanly
      problem = 'crash' -> a worker died (BrokenProcessPool)
      problem = 'stall' -> no task completed within `timeout`s => a worker is hung
    `suspects` = the cids still in-flight when the problem hit (bounded by window).
    """
    ex = ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=25)
    it = iter(list(remaining))
    inflight = {}  # fut -> cid
    problem = None

    def _submit_next():
        """Submit the next item. Returns False on iterator exhaustion OR a broken
        pool. ex.submit() itself raises BrokenProcessPool once a worker has died,
        so it MUST be guarded exactly like fut.result()/wait() — not doing so was
        the crash at the priming/refill submit sites."""
        nonlocal problem
        p = next(it, None)
        if p is None:
            return False
        try:
            inflight[ex.submit(work, p)] = p[0]
        except BrokenProcessPool:
            problem = "crash"
            return False
        return True

    try:
        for _ in range(workers):               # prime the window
            if not _submit_next():
                break
        while inflight and problem is None:
            done, _pending = wait(list(inflight), timeout=timeout,
                                  return_when=FIRST_COMPLETED)
            if not done:                       # nothing finished in `timeout`s -> hang
                problem = "stall"
                _kill_pool(ex)
                break
            for fut in done:
                cid = inflight.pop(fut)
                try:
                    rec = fut.result()
                except BrokenProcessPool:
                    problem = "crash"
                    inflight[fut] = cid   # retain the culprit as a suspect to isolate
                    break                 # (else a lone crashing task livelocks)
                except Exception as e:  # noqa: BLE001
                    rec = {"id": cid, "status": "timeout_or_crash", "error": str(e)[:160]}
                sink(rec)
                drop(cid)
                _submit_next()                 # refill; sets problem='crash' if pool broke
                if problem:
                    break
    finally:
        ex.shutdown(wait=False, cancel_futures=True)
    return list(inflight.values()), problem


def _isolate(suspects, payload_by_cid, timeout, sink, drop):
    """Re-run each suspect ALONE with a hard deadline; quarantine crash/hang culprits.
    workers=1 + max_tasks_per_child=1 => the one in-flight task IS the culprit, so a
    hang or crash is attributed exactly (no collateral)."""
    print(f"[resilient] isolating {len(suspects)} suspect(s) "
          f"(workers=1, deadline={timeout}s)", flush=True)
    for cid in suspects:
        p = payload_by_cid.get(cid)
        if p is None:
            continue
        ex = ProcessPoolExecutor(max_workers=1, max_tasks_per_child=1)
        fut = ex.submit(work, p)
        done, _ = wait([fut], timeout=timeout)
        if not done:                           # this exact complex hangs
            _kill_pool(ex)
            sink({"id": cid, "status": "crash_quarantine",
                  "error": f"isolated: hang > {timeout}s (native xtb/Architector loop); quarantined"})
        else:
            try:
                rec = fut.result()
            except BrokenProcessPool:
                rec = {"id": cid, "status": "crash_quarantine",
                       "error": "isolated: worker crash (segfault/abort); quarantined"}
            except Exception as e:  # noqa: BLE001
                rec = {"id": cid, "status": "timeout_or_crash", "error": str(e)[:160]}
            sink(rec)
        drop(cid)
        ex.shutdown(wait=False, cancel_futures=True)


def run_resilient(payloads, workers, timeout, sink):
    """Crash- AND hang-resilient pool runner.

    The first fix survived worker *crashes* (BrokenProcessPool) but not *hangs*:
    a complex whose native xtb/Architector build loops forever never completes,
    so the pool waits on it indefinitely (the 6,124 freeze). Here a fast
    sliding-window pass keeps the other workers busy; when it can make no progress
    within `timeout`s (stall) or a worker dies (crash), we tear the pool down
    (force-killing hung children) and isolate ONLY the in-flight suspects
    one-at-a-time with a hard deadline to pinpoint+quarantine the culprit. Repeats
    until every complex is built, failed, or quarantined.
    """
    payload_by_cid = {p[0]: p for p in payloads}
    state = {"remaining": list(payloads)}

    def drop(cid):
        state["remaining"] = [p for p in state["remaining"] if p[0] != cid]

    while state["remaining"]:
        before = len(state["remaining"])
        suspects, problem = _fast_pass(state["remaining"], workers, timeout, sink, drop)
        if problem is None:
            break
        print(f"[resilient] {problem}: {len(state['remaining'])} left, "
              f"{len(suspects)} suspect(s) to isolate", flush=True)
        if suspects:
            _isolate(suspects, payload_by_cid, timeout, sink, drop)
        elif len(state["remaining"]) == before:
            # crash/stall left no isolatable suspect AND nothing drained: the head
            # item is the culprit; force-quarantine it so we can't livelock on it.
            stuck = state["remaining"][0][0]
            print(f"[resilient] no-progress guard -> force-quarantine #{stuck}", flush=True)
            sink({"id": stuck, "status": "crash_quarantine",
                  "error": "force-quarantined: worker crashed with no isolatable suspect"})
            drop(stuck)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--limit", type=int, default=100)
    ap.add_argument("--ids", default="")
    ap.add_argument("--metal", default="")
    ap.add_argument("--all", action="store_true", help="all compounds, every metal")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=600)  # stall/hang deadline (tail tasks ≤~140s)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out, "struct"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "img"), exist_ok=True)

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    cols = "id, metal, oxidation_state, charge_complex, smiles_ligands, donor_atoms"
    if args.ids:
        ids = [int(x) for x in args.ids.split(",") if x.strip()]
        rows = con.execute(f"SELECT {cols} FROM complexes WHERE id IN ({','.join('?'*len(ids))})", ids).fetchall()
    elif args.all:
        rows = con.execute(f"SELECT {cols} FROM complexes ORDER BY id").fetchall()
    elif args.metal:
        rows = con.execute(f"SELECT {cols} FROM complexes WHERE metal=? ORDER BY id", (args.metal,)).fetchall()
    else:
        rows = con.execute(f"SELECT {cols} FROM complexes WHERE metal='Ir' ORDER BY id LIMIT ?",
                           (args.limit,)).fetchall()
    con.close()

    # resume: reload already-completed records from the JSONL checkpoint
    jsonl = os.path.join(args.out, "records.jsonl")
    done = {}
    if os.path.exists(jsonl):
        for line in open(jsonl):
            try:
                r = json.loads(line); done[r["id"]] = r
            except Exception:
                pass
    todo = [r for r in rows if r[0] not in done]
    print(f"[run] {len(rows)} complexes | {len(done)} done | {len(todo)} to build", flush=True)

    # precompute pydentate oracle for all unique ligand fragments (cached -> resumable)
    import coord_oracle as O
    frags = set()
    for r in rows:
        for f in r[4].split("."):
            f = f.strip()
            if f:
                frags.add(_canon(f))
    print(f"[run] precomputing oracle for {len(frags)} unique ligands...", flush=True)
    oracle = O.precompute(sorted(frags), out_path=os.path.join(args.out, "oracle.json"))

    payloads = [(r[0], r[1], r[2], r[3], r[4], r[5], oracle, args.out) for r in todo]
    jf = open(jsonl, "a")
    total = len(rows)
    counter = {"n": len(done)}

    def sink(rec):
        done[rec["id"]] = rec
        jf.write(json.dumps(rec) + "\n"); jf.flush()
        counter["n"] += 1
        print(f"[{counter['n']}/{total}] #{rec['id']} {rec.get('status'):<16} "
              f"isomers={len(rec.get('isomers', []))} {rec.get('build_s','')}s", flush=True)

    run_resilient(payloads, args.workers, args.timeout, sink)
    jf.close()

    import collections
    records = sorted(done.values(), key=lambda r: r["id"])
    manifest = {"n": len(records),
                "n_ok": sum(1 for r in records if r.get("status") == "ok"),
                "n_isomers_total": sum(len(r.get("isomers", [])) for r in records),
                "n_valid_struct": sum(1 for r in records for i in r.get("isomers", []) if i.get("valid")),
                "by_metal": dict(collections.Counter(r.get("metal") for r in records)),
                "records": records}
    json.dump(manifest, open(os.path.join(args.out, "manifest.json"), "w"), indent=2)
    print(f"[run] done. ok={manifest['n_ok']}/{manifest['n']} "
          f"structures={manifest['n_isomers_total']} valid={manifest['n_valid_struct']}", flush=True)


if __name__ == "__main__":
    main()
