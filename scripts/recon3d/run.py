"""Pilot/batch runner for 3D reconstruction.

Usage (inside arch-env):
  python run.py --db pilot.sqlite --sample 30 --mode fast --out out/pilot
  python run.py --db pilot.sqlite --ids 1,2,4820 --mode xtb --out out/xtb
Produces per-complex XYZ/MOL2/SDF + a QC report (JSON + Markdown).
"""
from __future__ import annotations
import argparse, json, os, sqlite3, sys, time, traceback
import multiprocessing as mp
from collections import Counter, defaultdict

import assign as A
import build as B
import validate as V


def _select(db, ids, sample, metal):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    cur = con.cursor()
    cols = "id, metal, oxidation_state, charge_complex, smiles_ligands, donor_atoms"
    if ids:
        q = f"SELECT {cols} FROM complexes WHERE id IN ({','.join('?'*len(ids))})"
        rows = cur.execute(q, ids).fetchall()
    elif metal:
        rows = cur.execute(f"SELECT {cols} FROM complexes WHERE metal=? ORDER BY id LIMIT ?",
                           (metal, sample)).fetchall()
    else:
        # stratified across metals, evenly spaced by id within each metal
        rows = []
        metals = [m for (m,) in cur.execute("SELECT DISTINCT metal FROM complexes")]
        per = max(1, sample // max(1, len(metals)))
        for m in metals:
            mr = cur.execute(f"SELECT {cols} FROM complexes WHERE metal=? ORDER BY id", (m,)).fetchall()
            if not mr:
                continue
            step = max(1, len(mr) // per)
            rows.extend(mr[::step][:per])
        rows = rows[:sample] if sample else rows
    con.close()
    return rows


def _work(row, mode, outdir, q):
    cid, metal, ox, charge, smiles, donor_json = row
    rec = {"id": cid, "metal": metal, "ox": ox, "charge_complex": charge,
           "smiles_ligands": smiles}
    try:
        donor_atoms = json.loads(donor_json) if donor_json else None
        asg = A.assign(metal, ox, smiles, donor_atoms)
        rec.update({
            "target_comp": dict(asg.target_comp), "target_cn": asg.target_cn,
            "pref_cn": asg.pref_cn, "geometry": asg.geometry, "spin": asg.spin,
            "n_ligand_units": len(asg.units), "assign_conf": asg.confidence,
            "assign_flags": asg.flags, "spectators": asg.spectators,
            "ligands": [{"smiles": u.smiles, "coordList": u.pocket,
                         "comp": dict(u.pocket_comp)} for u in asg.units],
        })
        t = time.time()
        res = B.build(asg, mode=mode)
        rec["build_status"] = res.get("status")
        rec["build_s"] = round(time.time() - t, 1)
        if res.get("status") == "ok":
            rec["method"] = res["method"]
            rec["n_atoms"] = res["n_atoms"]
            rec["built_charge"] = res.get("total_charge")
            passed, gates, neigh = V.validate(res, asg)
            rec["valid"] = passed
            rec["gates"] = gates
            paths = B.write_outputs(res, asg, outdir, cid)
            rec["files"] = paths
        else:
            rec["valid"] = False
            rec["build_error"] = res.get("error")
    except Exception as e:  # noqa: BLE001
        rec["build_status"] = "exception"
        rec["valid"] = False
        rec["exception"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    q.put(rec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--ids", default="")
    ap.add_argument("--sample", type=int, default=30)
    ap.add_argument("--metal", default="")
    ap.add_argument("--mode", default="fast", choices=["fast", "xtb"])
    ap.add_argument("--out", default="out/pilot")
    ap.add_argument("--timeout", type=int, default=0)
    args = ap.parse_args()

    ids = [int(x) for x in args.ids.split(",") if x.strip()] if args.ids else None
    rows = _select(args.db, ids, args.sample, args.metal or None)
    timeout = args.timeout or (600 if args.mode == "xtb" else 300)
    os.makedirs(args.out, exist_ok=True)

    print(f"[run] {len(rows)} complexes  mode={args.mode}  timeout={timeout}s  out={args.out}", flush=True)
    records = []
    ctx = mp.get_context("fork")
    for i, row in enumerate(rows, 1):
        q = ctx.Queue()
        p = ctx.Process(target=_work, args=(row, args.mode, args.out, q))
        p.start(); p.join(timeout)
        if p.is_alive():
            p.terminate(); p.join()
            rec = {"id": row[0], "metal": row[1], "ox": row[2],
                   "build_status": "timeout", "valid": False}
        else:
            rec = q.get() if not q.empty() else {"id": row[0], "metal": row[1],
                                                 "build_status": "died", "valid": False}
        records.append(rec)
        flag = "OK " if rec.get("valid") else "xx "
        print(f"[{i:>3}/{len(rows)}] {flag} #{rec['id']:<5} {rec.get('metal'):<3}"
              f"({rec.get('ox')}) {rec.get('build_status'):<10} "
              f"conf={rec.get('assign_conf','-'):<6} cn={rec.get('gates',{}).get('cn_got','-')}"
              f"/{rec.get('target_cn','-')} {rec.get('build_s','')}s "
              f"{rec.get('method','')}", flush=True)

    # ---- QC report ----
    rep = _report(records, args)
    with open(os.path.join(args.out, "qc.json"), "w") as f:
        json.dump({"args": vars(args), "summary": rep, "records": records}, f, indent=2)
    md = _report_md(rep, args)
    with open(os.path.join(args.out, "qc.md"), "w") as f:
        f.write(md)
    print("\n" + md, flush=True)


def _report(records, args):
    n = len(records)
    valid = [r for r in records if r.get("valid")]
    built = [r for r in records if r.get("build_status") == "ok"]
    by_metal = defaultdict(lambda: [0, 0])  # built_ok_valid, total
    for r in records:
        by_metal[r.get("metal")][1] += 1
        if r.get("valid"):
            by_metal[r.get("metal")][0] += 1
    conf = Counter(r.get("assign_conf") for r in records if r.get("assign_conf"))
    statuses = Counter(r.get("build_status") for r in records)
    comp_match = sum(1 for r in built if r.get("gates", {}).get("comp_match"))
    cn_match = sum(1 for r in built if r.get("gates", {}).get("cn_match"))
    flags = Counter(f.split(":")[0] for r in records for f in r.get("assign_flags", []))
    times = [r["build_s"] for r in built if "build_s" in r]
    return {
        "n": n, "n_built": len(built), "n_valid": len(valid),
        "valid_rate": round(len(valid) / n, 3) if n else 0,
        "build_rate": round(len(built) / n, 3) if n else 0,
        "by_metal": {k: {"valid": v[0], "total": v[1]} for k, v in by_metal.items()},
        "assign_confidence": dict(conf), "build_status": dict(statuses),
        "comp_match_of_built": f"{comp_match}/{len(built)}",
        "cn_match_of_built": f"{cn_match}/{len(built)}",
        "assign_flags": dict(flags),
        "build_time_s": {"mean": round(sum(times)/len(times), 1) if times else 0,
                         "max": max(times) if times else 0},
    }


def _report_md(rep, args):
    L = [f"# 3D reconstruction QC — mode={args.mode}, n={rep['n']}", ""]
    L.append(f"- built (geometry produced): **{rep['n_built']}/{rep['n']}** ({rep['build_rate']:.0%})")
    L.append(f"- valid (passed all gates):  **{rep['n_valid']}/{rep['n']}** ({rep['valid_rate']:.0%})")
    L.append(f"- donor-composition match (of built): {rep['comp_match_of_built']}")
    L.append(f"- CN match (of built): {rep['cn_match_of_built']}")
    L.append(f"- build time: mean {rep['build_time_s']['mean']}s, max {rep['build_time_s']['max']}s")
    L.append("")
    L.append("| metal | valid/total |"); L.append("|---|---|")
    for m, v in sorted(rep["by_metal"].items()):
        L.append(f"| {m} | {v['valid']}/{v['total']} |")
    L.append("")
    L.append(f"- assign confidence: {rep['assign_confidence']}")
    L.append(f"- build status: {rep['build_status']}")
    L.append(f"- assign flags: {rep['assign_flags']}")
    return "\n".join(L)


if __name__ == "__main__":
    main()
