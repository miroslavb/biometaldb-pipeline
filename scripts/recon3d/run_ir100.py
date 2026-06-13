"""Reconstruct the first N iridium complexes -> stereoisomers + TREXIO + manifest.

Pipeline per complex: pydentate oracle -> assign -> build_isomers (cis/trans,
fac/mer) -> validate each -> write XYZ/mol2/TREXIO + 2D PNG -> manifest record.
The manifest feeds review_ui.py (manual-verification interface).

Run (on hive or locally, inside arch-env):
  python run_ir100.py --db pilot.sqlite --limit 100 --workers 16 --out out/ir100
"""
from __future__ import annotations
import argparse, json, os, sqlite3, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--limit", type=int, default=100)
    ap.add_argument("--ids", default="")
    ap.add_argument("--metal", default="")
    ap.add_argument("--all", action="store_true", help="all compounds, every metal")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=900)
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
    n = len(done)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(work, p): p[0] for p in payloads}
        for fut in as_completed(futs):
            cid = futs[fut]
            try:
                rec = fut.result(timeout=args.timeout)
            except Exception as e:  # noqa: BLE001
                rec = {"id": cid, "status": "timeout_or_crash", "error": str(e)[:120]}
            done[cid] = rec
            jf.write(json.dumps(rec) + "\n"); jf.flush()
            n += 1
            print(f"[{n}/{len(rows)}] #{rec['id']} {rec.get('status'):<12} "
                  f"isomers={len(rec.get('isomers', []))} {rec.get('build_s','')}s", flush=True)
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
