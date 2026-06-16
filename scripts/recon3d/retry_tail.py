"""Targeted re-run of the FAILED tail of the BiometalDB 3D reconstruction.

The full run (run_ir100.py --all) left 353/9414 complexes without a structure:
~330 `no_structure` (Architector produced no conformer) + ~23 `crash_quarantine`
(native xtb/Architector hang or segfault). This module re-attacks ONLY those,
with an escalating strategy ladder that goes beyond the original build ladder,
and — for the ones that still cannot be built — attaches a clear, human-readable
failure reason that the review UI surfaces verbatim.

Strategy ladder (first success wins; each is a genuinely different attempt):
  S1  xtb  + donor coordList + 12 conformers      (stochastic re-roll of the embed)
  S2  UFF  + donor coordList + 12 conformers      (permissive placement, no xtb)
  S3  xtb  + Architector auto-donor (no coordList) (our pocket may be wrong)
  S4  UFF  + Architector auto-donor (no coordList) (last resort: any sane geometry)

What it does NOT do: invent chemistry. If every strategy fails the record stays
`no_structure` and gets fail_reason/fail_label/fail_hint so the UI states WHY.

Run on hive inside arch-env:
  python retry_tail.py --db pilot.sqlite --out out/full --workers 24 --timeout 900
Resumable: reloads out/full/retry.jsonl; merges into records.jsonl + manifest.json.
"""
from __future__ import annotations
import argparse, ast, collections, json, os, re, time, traceback
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from concurrent.futures.process import BrokenProcessPool

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

import assign as A
import build as B
import validate as V
import trexio_writer as TX
from run_ir100 import render_png, _canon


# ---------------------------------------------------------------------------
# build strategies
# ---------------------------------------------------------------------------
def _dedup_build(asg, relax, n_conf, use_coordlist, max_isomers=2):
    """One build attempt -> list of distinct-stereoisomer result dicts (or [])."""
    from architector import build_complex
    if not asg.units:
        return []
    amethod = "GFN2-xTB" if relax else "UFF"

    def ligs():
        out = []
        for u in asg.units:
            d = {"smiles": u.smiles}
            if use_coordlist and u.pocket:
                d["coordList"] = list(u.pocket)
                if u.ligType:
                    d["ligType"] = u.ligType
            out.append(d)
        return out

    params = dict(metal_ox=asg.ox, metal_spin=asg.spin, return_only_1=False,
                  n_conformers=n_conf, n_symmetries=max(n_conf + 6, 14),
                  save_init_geos=False, debug=False, relax=relax,
                  assemble_method=amethod)
    inp = {"core": {"metal": asg.metal, "coreType": asg.geometry},
           "ligands": ligs(), "parameters": params}
    try:
        out = build_complex(inp)
    except Exception:
        return []
    if not out:
        return []
    confs = [(k, v) for k, v in out.items()
             if "_init_only" not in k and v.get("ase_atoms") is not None]
    confs.sort(key=lambda kv: kv[1].get("energy", 1e18))
    kept, seen = [], set()
    for _, v in confs:
        sig = B._stereo_signature(v["ase_atoms"], asg.metal)
        if sig in seen:
            continue
        seen.add(sig)
        r = B._mk_result(v, amethod + ("+relax" if relax else ""), 0)
        r["stereo_sig"] = str(sig)
        kept.append(r)
        if len(kept) >= max_isomers:
            break
    for i, r in enumerate(kept):
        r["isomer"] = "only" if len(kept) == 1 else f"isomer{i + 1}"
    return kept


# (name, relax, n_conf, use_coordlist)
_LADDER = [
    ("S1_xtb_donor_12conf", True,  12, True),
    ("S2_uff_donor_12conf", False, 12, True),
    ("S3_xtb_auto_donor",   True,  12, False),
    ("S4_uff_auto_donor",   False, 12, False),
]


def _try_strategies(asg):
    """Run the ladder; return (isomers, strategy_name, tried_list)."""
    tried = []
    for name, relax, n_conf, use_cl in _LADDER:
        tried.append(name)
        isos = _dedup_build(asg, relax=relax, n_conf=n_conf, use_coordlist=use_cl)
        if isos:
            return isos, name, tried
    return [], None, tried


# ---------------------------------------------------------------------------
# failure classification (only reached when the whole ladder fails)
# ---------------------------------------------------------------------------
_REASONS = {
    "native_runtime_failure": (
        "Builder crash / infinite loop",
        "The xtb/Architector native build crashes or hangs even when run in "
        "isolation. Not recoverable without a code-level fix to the builder."),
    "corrupt_source_donors": (
        "Implausible donor data in source DB",
        "The donor-atom record from the source database is chemically "
        "implausible (e.g. dozens of O/S donors on one ligand), so no real "
        "coordination sphere can be assembled. The input SMILES / donor entry "
        "must be corrected at the source."),
    "unparseable_smiles": (
        "Unparseable ligand SMILES",
        "At least one ligand SMILES could not be parsed by RDKit, so the "
        "ligand could not be placed. A corrected structure is required."),
    "insufficient_donors": (
        "Coordination sphere under-filled",
        "Too few donor atoms could be identified to fill the metal's "
        "coordination sphere. Needs manual donor-atom curation."),
    "haptic_overcrowded": (
        "Haptic ligand has no room",
        "A haptic (Cp/arene) ligand cannot fit alongside the other ligands in "
        "the coordination sphere. The ligand set is likely over-specified."),
    "geometry_unassemblable": (
        "Geometry could not be assembled",
        "Architector could not pack the assigned donor set into any tried "
        "polyhedron (xtb and UFF, with and without donor hints, 12 conformers "
        "each). This usually means the donor/denticity assignment is wrong and "
        "needs manual review."),
}


def _has_absurd_want(flags):
    """A `comp_got{...}_want{...}` flag whose 'want' is chemically impossible."""
    for f in flags or []:
        m = re.search(r"_want(\{.*\})\s*$", str(f))
        if not m:
            continue
        try:
            want = ast.literal_eval(m.group(1))
        except Exception:
            continue
        if isinstance(want, dict) and want:
            if any(int(v) > 8 for v in want.values()) or sum(int(v) for v in want.values()) > 12:
                return True
    return False


def classify_failure(asg, smiles, donor_atoms):
    flags = list(asg.flags or [])
    has_unparsed = any(getattr(u, "role", None) == "unparsed" for u in asg.units) \
        or any("unparsed" in str(f) for f in flags)
    if has_unparsed:
        code = "unparseable_smiles"
    elif _has_absurd_want(flags):
        code = "corrupt_source_donors"
    elif any(str(f).startswith("haptic_no_room") for f in flags):
        code = "haptic_overcrowded"
    elif any(str(f).startswith("cn_unfilled") for f in flags) \
            or "no_donor_atoms_prior" in flags:
        code = "insufficient_donors"
    else:
        code = "geometry_unassemblable"
    label, hint = _REASONS[code]
    detail_flags = [f for f in flags if "comp_got" in str(f) or "donor_prior" in str(f)
                    or "cn_unfilled" in str(f) or "unused" in str(f)
                    or "haptic" in str(f)][:3]
    return {"fail_reason": code, "fail_label": label, "fail_hint": hint,
            "fail_flags": detail_flags}


# ---------------------------------------------------------------------------
# per-complex worker (mirrors run_ir100.work but with the retry ladder)
# ---------------------------------------------------------------------------
def work_retry(args):
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
        png_path = os.path.join(outdir, "img", f"c{cid}.png")
        rec["png"] = (os.path.basename(png_path) if os.path.exists(png_path)
                      else render_png(cid, metal, ox, smiles, png_path))

        t = time.time()
        isos, strat, tried = _try_strategies(asg)
        rec["build_s"] = round(time.time() - t, 1)

        if isos:
            for k, res in enumerate(isos):
                label = res.get("isomer", f"iso{k + 1}")
                passed, gates, neigh = V.validate(res, asg)
                stem = os.path.join(outdir, "struct", f"c{cid}_{label}")
                paths = B.write_outputs(res, asg, os.path.join(outdir, "struct"), f"{cid}_{label}")
                tx = None
                try:
                    tx = TX.write_trexio(stem + ".h5", res["symbols"], res["positions"],
                                         res.get("total_charge"), res.get("n_unpaired"),
                                         title=f"BiometalDB #{cid} {metal}({ox}) {label}")
                except Exception as e:  # noqa: BLE001
                    rec.setdefault("trexio_errors", []).append(
                        f"{label}: {type(e).__name__} {str(e)[:80]}")
                rec["isomers"].append({
                    "label": label, "valid": passed, "energy": res.get("energy"),
                    "total_charge": res.get("total_charge"),
                    "mult": (res.get("n_unpaired") or 0) + 1,
                    "n_atoms": res.get("n_atoms"), "method": res.get("method"),
                    "stereo_sig": res.get("stereo_sig"),
                    "metal_neighbors": [[s, d] for _, s, d in neigh],
                    "gates": {k2: gates[k2] for k2 in ("no_clash", "bond_lengths_ok",
                              "min_heavy_dist", "cn_got", "comp_got") if k2 in gates},
                    "files": {**{k2: os.path.basename(v) for k2, v in paths.items()},
                              **({"trexio": os.path.basename(tx)} if tx else {})},
                })
            rec["status"] = "ok"
            rec["retry_meta"] = {"recovered_by": strat, "strategies_tried": tried}
        else:
            rec["status"] = "no_structure"
            rec.update(classify_failure(asg, smiles, donor_atoms))
            rec["retry_meta"] = {"recovered_by": None, "strategies_tried": tried}
    except Exception as e:  # noqa: BLE001
        rec["status"] = "no_structure"
        rec["fail_reason"] = "geometry_unassemblable"
        rec["fail_label"] = _REASONS["geometry_unassemblable"][0]
        rec["fail_hint"] = (_REASONS["geometry_unassemblable"][1]
                            + f"  (assign/build raised {type(e).__name__})")
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-500:]
        rec["retry_meta"] = {"recovered_by": None, "strategies_tried": []}
    return rec


# ---------------------------------------------------------------------------
# resilient pool (crash + hang isolation) — local copy bound to work_retry
# ---------------------------------------------------------------------------
def _kill_pool(ex):
    for proc in list(getattr(ex, "_processes", {}).values()):
        try:
            proc.kill()
        except Exception:  # noqa: BLE001
            pass


def _payload_stub(p, verdict):
    """Quarantine record carrying the payload's metadata (metal/ox/smiles +
    depiction if the worker already rendered one) so a crash/hang card is still
    informative — not a blank record. `p` is the payload tuple or None."""
    rec = {"isomers": []}
    if p:
        cid, metal, ox, charge, smiles = p[0], p[1], p[2], p[3], p[4]
        outdir = p[7]
        rec.update({"id": cid, "metal": metal, "ox": ox, "charge_complex": charge,
                    "smiles_ligands": smiles})
        png = os.path.join(outdir, "img", f"c{cid}.png")
        if os.path.exists(png):
            rec["png"] = os.path.basename(png)
    rec.update(verdict)
    rec.setdefault("retry_meta", {"recovered_by": None, "strategies_tried": ["isolated"]})
    return rec


def _native_verdict(suffix=""):
    return {"status": "no_structure", "fail_reason": "native_runtime_failure",
            "fail_label": _REASONS["native_runtime_failure"][0],
            "fail_hint": _REASONS["native_runtime_failure"][1] + suffix}


def _fast_pass(remaining, workers, timeout, sink, drop):
    ex = ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=20)
    pbc = {p[0]: p for p in remaining}
    it = iter(list(remaining))
    inflight, problem = {}, None

    def _submit_next():
        nonlocal problem
        p = next(it, None)
        if p is None:
            return False
        try:
            inflight[ex.submit(work_retry, p)] = p[0]
        except BrokenProcessPool:
            problem = "crash"
            return False
        return True

    try:
        for _ in range(workers):
            if not _submit_next():
                break
        while inflight and problem is None:
            done, _pending = wait(list(inflight), timeout=timeout, return_when=FIRST_COMPLETED)
            if not done:
                problem = "stall"
                _kill_pool(ex)
                break
            for fut in done:
                cid = inflight.pop(fut)
                try:
                    rec = fut.result()
                except BrokenProcessPool:
                    problem = "crash"
                    inflight[fut] = cid   # keep the culprit in the suspect set
                    break                 # (else a lone crashing task livelocks)
                except Exception as e:  # noqa: BLE001
                    rec = _payload_stub(pbc.get(cid), {**_native_verdict(), "error": str(e)[:160]})
                sink(rec)
                drop(cid)
                _submit_next()
                if problem:
                    break
    finally:
        ex.shutdown(wait=False, cancel_futures=True)
    return list(inflight.values()), problem


def _isolate(suspects, payload_by_cid, timeout, sink, drop):
    print(f"[retry] isolating {len(suspects)} suspect(s) (workers=1, deadline={timeout}s)",
          flush=True)
    for cid in suspects:
        p = payload_by_cid.get(cid)
        if p is None:
            continue
        ex = ProcessPoolExecutor(max_workers=1, max_tasks_per_child=1)
        fut = ex.submit(work_retry, p)
        done, _ = wait([fut], timeout=timeout)
        if not done:
            _kill_pool(ex)
            sink(_payload_stub(p, _native_verdict(f"  (isolated hang > {timeout}s)")))
        else:
            try:
                rec = fut.result()
            except BrokenProcessPool:
                rec = _payload_stub(p, _native_verdict("  (isolated worker segfault/abort)"))
            except Exception as e:  # noqa: BLE001
                rec = _payload_stub(p, {**_native_verdict(), "error": str(e)[:160]})
            sink(rec)
        drop(cid)
        ex.shutdown(wait=False, cancel_futures=True)


def run_resilient(payloads, workers, timeout, sink):
    payload_by_cid = {p[0]: p for p in payloads}
    state = {"remaining": list(payloads)}

    def drop(cid):
        state["remaining"] = [p for p in state["remaining"] if p[0] != cid]

    while state["remaining"]:
        before = len(state["remaining"])
        suspects, problem = _fast_pass(state["remaining"], workers, timeout, sink, drop)
        if problem is None:
            break
        print(f"[retry] {problem}: {len(state['remaining'])} left, "
              f"{len(suspects)} suspect(s) to isolate", flush=True)
        if suspects:
            _isolate(suspects, payload_by_cid, timeout, sink, drop)
        elif len(state["remaining"]) == before:
            # crash/stall left no isolatable suspect AND nothing drained: the head
            # item is the culprit (e.g. a native abort on submit). Force-quarantine
            # it so the loop can never spin forever on one bad complex.
            stuck_p = state["remaining"][0]
            stuck = stuck_p[0]
            print(f"[retry] no-progress guard -> force-quarantine #{stuck}", flush=True)
            sink(_payload_stub(stuck_p, _native_verdict(
                "  (force-quarantined: worker crashed with no isolatable suspect)")))
            drop(stuck)


# ---------------------------------------------------------------------------
# merge retry results INTO the (enriched) manifest.json + records.jsonl
# ---------------------------------------------------------------------------
def _normalize_crash(r):
    if r.get("status") in ("crash_quarantine", "timeout_or_crash", "exception") \
            and not r.get("fail_reason"):
        r["status"] = "no_structure"
        r["fail_reason"] = "native_runtime_failure"
        r["fail_label"] = _REASONS["native_runtime_failure"][0]
        r["fail_hint"] = _REASONS["native_runtime_failure"][1]
    return r


def _recompute_manifest_counts(man):
    recs = man["records"]
    man["n"] = len(recs)
    man["n_ok"] = sum(1 for r in recs if r.get("status") == "ok")
    man["n_isomers_total"] = sum(len(r.get("isomers", [])) for r in recs)
    man["n_valid_struct"] = sum(1 for r in recs for i in r.get("isomers", []) if i.get("valid"))
    man["n_recovered"] = sum(1 for r in recs if (r.get("retry_meta") or {}).get("recovered_by"))
    man["n_still_failed"] = sum(1 for r in recs if r.get("status") != "ok")
    man["by_metal"] = dict(collections.Counter(r.get("metal") for r in recs))
    man["by_fail_reason"] = dict(collections.Counter(
        r.get("fail_reason") for r in recs if r.get("status") != "ok"))


def _enrich_recovered(out_dir, recovered_ids):
    """Best-effort: give recovered records the same enrichment as the rest
    (enantiomer mirror, TREXIO text, per-compound zip). The post-passes are
    idempotent (skip already-enriched records), so running them over the full
    manifest only touches the freshly recovered ones. Never raises."""
    if not recovered_ids:
        return "skipped (none recovered)"
    import subprocess
    here = os.path.dirname(os.path.abspath(__file__))
    logs = []
    for script in ("add_enantiomers.py", "add_trexio_text.py", "make_archives.py"):
        try:
            p = subprocess.run([__import__("sys").executable,
                                os.path.join(here, script), out_dir],
                               capture_output=True, text=True, timeout=3600)
            logs.append(f"{script}: rc={p.returncode} {p.stdout.strip()[-120:]}")
        except Exception as e:  # noqa: BLE001
            logs.append(f"{script}: ERROR {type(e).__name__} {str(e)[:80]}")
    return " | ".join(logs)


def merge_and_manifest(out_dir):
    manifest_path = os.path.join(out_dir, "manifest.json")
    records_path = os.path.join(out_dir, "records.jsonl")
    retry_path = os.path.join(out_dir, "retry.jsonl")

    retry_recs = {}
    for line in open(retry_path):
        line = line.strip()
        if not line:
            continue
        try:
            r = _normalize_crash(json.loads(line))
        except Exception:
            continue
        retry_recs[r["id"]] = r

    # 1) merge into the ENRICHED manifest (preserve enantiomers/archives/trexio_txt
    #    for the 9061 already-built complexes; replace only the retried ids)
    man = json.load(open(manifest_path))
    idx = {r["id"]: i for i, r in enumerate(man["records"])}
    n_over = 0
    for cid, r in retry_recs.items():
        is_stub = ("smiles_ligands" not in r) and ("ligands" not in r)
        if cid in idx and is_stub:
            # bare isolation stub (crash/hang): keep the original card's metadata
            # (smiles, depiction, ligand table, geometry) and overlay only the verdict
            merged = dict(man["records"][idx[cid]])
            for k in ("status", "fail_reason", "fail_label", "fail_hint",
                      "fail_flags", "retry_meta", "error", "trace"):
                if k in r:
                    merged[k] = r[k]
            merged["isomers"] = []
            man["records"][idx[cid]] = merged
        elif cid in idx:
            man["records"][idx[cid]] = r
        else:
            man["records"].append(r)
            idx[cid] = len(man["records"]) - 1
        n_over += 1
    man["records"] = [_normalize_crash(r) for r in man["records"]]
    man["records"].sort(key=lambda r: r["id"])

    bak = manifest_path + ".bak.pre-retry"
    if not os.path.exists(bak):
        import shutil
        shutil.copy2(manifest_path, bak)
    _recompute_manifest_counts(man)
    json.dump(man, open(manifest_path, "w"), indent=2)

    # 2) enrich the freshly recovered records (idempotent passes), then refresh counts
    recovered_ids = {r["id"] for r in retry_recs.values()
                     if (r.get("retry_meta") or {}).get("recovered_by")}
    enrich_log = _enrich_recovered(out_dir, recovered_ids)
    man = json.load(open(manifest_path))
    _recompute_manifest_counts(man)
    json.dump(man, open(manifest_path, "w"), indent=2)

    # 3) keep records.jsonl (the resume ledger) consistent — replace retried ids
    if os.path.exists(records_path):
        by_id = {}
        for line in open(records_path):
            line = line.strip()
            if not line:
                continue
            try:
                rr = json.loads(line)
                by_id[rr["id"]] = rr
            except Exception:
                pass
        for cid, r in retry_recs.items():
            by_id[cid] = r
        rbak = records_path + ".bak.pre-retry"
        if not os.path.exists(rbak):
            os.replace(records_path, rbak)
        with open(records_path, "w") as f:
            for rr in sorted(by_id.values(), key=lambda r: r["id"]):
                f.write(json.dumps(rr) + "\n")

    man["_enrich_log"] = enrich_log
    return man, n_over


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", default="out/full")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--ids", default="", help="comma-list to restrict (debug)")
    ap.add_argument("--merge-only", action="store_true",
                    help="skip compute, just merge existing retry.jsonl + manifest")
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out, "struct"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "img"), exist_ok=True)

    if args.merge_only:
        man, n = merge_and_manifest(args.out)
        print(f"[retry] merge-only: applied {n} retry records | ok={man['n_ok']}/{man['n']} "
              f"recovered={man['n_recovered']} still_failed={man['n_still_failed']}", flush=True)
        return

    # failed ids from the current records.jsonl (no isomers)
    records_path = os.path.join(args.out, "records.jsonl")
    failed_ids = []
    for line in open(records_path):
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except Exception:
            continue
        if not (r.get("isomers") or []):
            failed_ids.append(r["id"])
    if args.ids:
        want = {int(x) for x in args.ids.split(",") if x.strip()}
        failed_ids = [i for i in failed_ids if i in want]
    failed_ids = sorted(set(failed_ids))

    # already-retried ids (resume)
    retry_path = os.path.join(args.out, "retry.jsonl")
    done = set()
    if os.path.exists(retry_path):
        for line in open(retry_path):
            try:
                done.add(json.loads(line)["id"])
            except Exception:
                pass
    todo_ids = [i for i in failed_ids if i not in done]
    print(f"[retry] {len(failed_ids)} failed | {len(done)} already retried | "
          f"{len(todo_ids)} to attempt", flush=True)

    if todo_ids:
        import sqlite3
        con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
        cols = "id, metal, oxidation_state, charge_complex, smiles_ligands, donor_atoms"
        rows = con.execute(
            f"SELECT {cols} FROM complexes WHERE id IN ({','.join('?' * len(todo_ids))})",
            todo_ids).fetchall()
        con.close()

        # reuse the cached pydentate oracle (don't recompute)
        oracle_path = os.path.join(args.out, "oracle.json")
        oracle = json.load(open(oracle_path)) if os.path.exists(oracle_path) else {}
        print(f"[retry] loaded oracle with {len(oracle)} cached ligands", flush=True)

        payloads = [(r[0], r[1], r[2], r[3], r[4], r[5], oracle, args.out) for r in rows]
        jf = open(retry_path, "a")
        counter = {"n": len(done)}
        total = len(failed_ids)

        def sink(rec):
            jf.write(json.dumps(rec) + "\n")
            jf.flush()
            counter["n"] += 1
            rb = (rec.get("retry_meta") or {}).get("recovered_by")
            tag = (f"RECOVERED:{rb}" if rb else f"FAIL:{rec.get('fail_reason', '?')}")
            print(f"[{counter['n']}/{total}] #{rec['id']} {rec.get('status'):<13} "
                  f"{tag} {rec.get('build_s', '')}s", flush=True)

        run_resilient(payloads, args.workers, args.timeout, sink)
        jf.close()

    man, n_over = merge_and_manifest(args.out)
    print(f"[retry] DONE. applied {n_over} retry records | ok={man['n_ok']}/{man['n']} "
          f"recovered={man['n_recovered']} still_failed={man['n_still_failed']}", flush=True)
    print(f"[retry] by_fail_reason: {man['by_fail_reason']}", flush=True)


if __name__ == "__main__":
    main()
