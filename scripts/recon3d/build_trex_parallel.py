"""build_trex_parallel — parallel T-REX generator using ProcessPoolExecutor.

Designed for hive3070t06 (56 cores, 128 GB RAM). Processes all 9414 complexes
in parallel batches, with per-complex timeout protection via subprocess.

Usage:
    cd /root/biometaldb-3d
    export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
    ./bin/micromamba run -p ./arch-env python3 polynuclear/build_trex_parallel.py
"""
import os
import sys
import json
import time
import tempfile
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter

D_BLOCK = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
           "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
           "Hf","Ta","W","Re","Os","Ir","Pt","Au"}

TREX_WORKER = """import sys, json
sys.path.insert(0, '/tmp/trex_repo')
from trex.xyz_to_trex import xyz_to_trex
args = json.loads(sys.argv[1])
result = xyz_to_trex(args['xyz_path'], overall_charge=args['charge'])
print(json.dumps({'ok': True, 'trex': result}))
"""


def normalize_xyz_for_trex(src_path, charge):
    """Create a temporary XYZ file with a T-REX-friendly comment line."""
    with open(src_path) as f:
        lines = f.readlines()
    n_atoms = lines[0].strip()
    new_comment = f"charge={charge}"
    fd, dst = tempfile.mkstemp(suffix=".xyz", prefix="trex_")
    with os.fdopen(fd, "w") as f:
        f.write(f"{n_atoms}\n{new_comment}\n")
        f.writelines(lines[2:])
    return dst


def count_metals_in_xyz(path):
    """Return number of d-block metal atoms in the XYZ file."""
    n_total = None
    n = 0
    with open(path) as f:
        for j, ln in enumerate(f):
            if j == 0:
                try:
                    n_total = int(ln.strip())
                except ValueError:
                    n_total = -1
                continue
            if j == 1:
                continue
            if not ln.strip():
                continue
            sym = ln.split()[0]
            if sym in D_BLOCK:
                n += 1
            if n_total > 0 and j - 1 >= n_total:
                break
    return n


def process_one(args_dict):
    """Process a single complex. Returns (cid, result_dict) or (cid, None) on failure.

    Runs in a worker process. Each call to xyz_to_trex is wrapped in a
    subprocess with 15s timeout to prevent hangs.
    """
    cid = args_dict["cid"]
    xyz_path = args_dict["xyz_path"]
    charge = args_dict["charge"]
    method = args_dict["method"]

    # Check polynuclear
    n_metals = count_metals_in_xyz(xyz_path)
    if n_metals != 1:
        return (cid, None, "polynuclear")

    tmp_xyz = None
    try:
        tmp_xyz = normalize_xyz_for_trex(xyz_path, charge)
        # Write worker script once per process
        worker_script = os.path.join(tempfile.gettempdir(), f"_trex_worker_{os.getpid()}.py")
        if not os.path.exists(worker_script):
            with open(worker_script, "w") as f:
                f.write(TREX_WORKER)

        proc = subprocess.run(
            [sys.executable, worker_script, json.dumps({"xyz_path": tmp_xyz, "charge": charge})],
            capture_output=True, text=True, timeout=15,
            env={**os.environ, "PYTHONPATH": "/tmp/trex_repo:" + os.environ.get("PYTHONPATH", "")}
        )
        if proc.returncode != 0:
            return (cid, None, f"exit={proc.returncode}: {proc.stderr[:100]}")
        result = json.loads(proc.stdout.strip())
        return (cid, {"trex": result["trex"], "charge": charge, "method": method}, "ok")
    except subprocess.TimeoutExpired:
        return (cid, None, "timeout")
    except Exception as e:
        return (cid, None, f"{type(e).__name__}: {str(e)[:80]}")
    finally:
        if tmp_xyz and os.path.exists(tmp_xyz):
            try:
                os.unlink(tmp_xyz)
            except:
                pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/root/biometaldb-3d/out/full")
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--batch-size", type=int, default=200)
    ap.add_argument("--skip-existing", action="store_true", default=True)
    args = ap.parse_args()

    manifest_path = os.path.join(args.out, "manifest.json")
    man = json.load(open(manifest_path))
    recs = man["records"]

    out_path = os.path.join(args.out, "trex.json")
    trex_data = {}
    if args.skip_existing and os.path.exists(out_path):
        trex_data = json.load(open(out_path))
        print(f"Loaded {len(trex_data)} existing T-REX entries", flush=True)

    # Build task list
    tasks = []
    n_skip_poly = 0
    n_skip_nostruct = 0
    for r in recs:
        cid = r["id"]
        if str(cid) in trex_data:
            continue
        if r.get("status") != "ok":
            n_skip_nostruct += 1
            continue
        isomers = r.get("isomers", [])
        iso = next((i for i in isomers if i.get("files", {}).get("xyz")), None)
        if not iso:
            n_skip_nostruct += 1
            continue
        xyz_fn = iso["files"]["xyz"]
        xyz_path = os.path.join(args.out, "struct", xyz_fn)
        if not os.path.exists(xyz_path):
            n_skip_nostruct += 1
            continue
        tasks.append({
            "cid": cid,
            "xyz_path": xyz_path,
            "charge": iso.get("total_charge", 0) or 0,
            "method": iso.get("method", ""),
        })

    print(f"Tasks to process: {len(tasks)} (skip_poly will be detected at runtime, "
          f"pre-skip: no_struct={n_skip_nostruct})", flush=True)

    n_ok = len(trex_data)
    n_fail = 0
    n_poly = 0
    n_timeout = 0
    t0 = time.time()

    # Process in batches
    for batch_start in range(0, len(tasks), args.batch_size):
        batch = tasks[batch_start:batch_start + args.batch_size]
        batch_end = min(batch_start + len(batch), len(tasks))

        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(process_one, t): t for t in batch}
            for fut in as_completed(futures):
                cid, result, status = fut.result()
                if status == "ok":
                    trex_data[str(cid)] = result
                    n_ok += 1
                elif status == "polynuclear":
                    n_poly += 1
                elif status == "timeout":
                    n_timeout += 1
                    n_fail += 1
                else:
                    n_fail += 1

        elapsed = time.time() - t0
        # Save progress
        with open(out_path, "w") as f:
            json.dump(trex_data, f, separators=(",", ":"))
        print(f"  {batch_end}/{len(tasks)}: ok={n_ok}, poly_skip={n_poly}, "
              f"fail={n_fail}, timeout={n_timeout} ({elapsed:.0f}s)",
              flush=True)

    elapsed = time.time() - t0
    print(f"DONE: {n_ok} T-REX generated, {n_poly} polynuclear skipped, "
          f"{n_fail} failed ({n_timeout} timeout) ({elapsed:.0f}s)")


if __name__ == "__main__":
    main()
