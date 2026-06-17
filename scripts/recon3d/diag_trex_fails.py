"""Diagnose T-REX failures: classify each missing ID by error type.

Usage on hive:
    cd /root/biometaldb-3d
    export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
    ./bin/micromamba run -p ./arch-env python3 polynuclear/diag_trex_fails.py
"""
import os, sys, json, time, tempfile, subprocess
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

def normalize_xyz(src_path, charge):
    with open(src_path) as f:
        lines = f.readlines()
    n_atoms = lines[0].strip()
    fd, dst = tempfile.mkstemp(suffix=".xyz", prefix="trex_")
    with os.fdopen(fd, "w") as f:
        f.write(f"{n_atoms}\ncharge={charge}\n")
        f.writelines(lines[2:])
    return dst

def count_metals(path):
    n_total = None
    n = 0
    with open(path) as f:
        for j, ln in enumerate(f):
            if j == 0:
                try: n_total = int(ln.strip())
                except: n_total = -1
                continue
            if j == 1: continue
            if not ln.strip(): continue
            if ln.split()[0] in D_BLOCK: n += 1
            if n_total > 0 and j - 1 >= n_total: break
    return n

def diag_one(args_dict):
    cid = args_dict["cid"]
    xyz_path = args_dict["xyz_path"]
    charge = args_dict["charge"]
    
    n_metals = count_metals(xyz_path)
    if n_metals != 1:
        return (cid, "polynuclear", f"{n_metals} metals", 0)
    
    # Get atom count
    with open(xyz_path) as f:
        n_atoms = int(f.readline().strip())
    
    tmp_xyz = None
    t0 = time.time()
    try:
        tmp_xyz = normalize_xyz(xyz_path, charge)
        worker_script = os.path.join(tempfile.gettempdir(), f"_tw_{os.getpid()}.py")
        if not os.path.exists(worker_script):
            with open(worker_script, "w") as f:
                f.write(TREX_WORKER)
        
        proc = subprocess.run(
            [sys.executable, worker_script, json.dumps({"xyz_path": tmp_xyz, "charge": charge})],
            capture_output=True, text=True, timeout=30,
            env={**os.environ, "PYTHONPATH": "/tmp/trex_repo:" + os.environ.get("PYTHONPATH", "")}
        )
        elapsed = time.time() - t0
        if proc.returncode != 0:
            stderr = proc.stderr[:200]
            if "Kekulize" in stderr:
                return (cid, "kekulize", stderr.strip(), elapsed)
            elif "Distance" in stderr or "Warning" in stderr:
                # Warnings only, check if stdout has result
                if proc.stdout.strip():
                    try:
                        result = json.loads(proc.stdout.strip())
                        return (cid, "ok_with_warnings", f"warnings, exit={proc.returncode}", elapsed)
                    except:
                        pass
                return (cid, "exit_with_warnings", stderr.strip()[:150], elapsed)
            else:
                return (cid, "other_exit", f"exit={proc.returncode}: {stderr[:150]}", elapsed)
        result = json.loads(proc.stdout.strip())
        return (cid, "ok", result["trex"][:60], elapsed)
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        return (cid, "timeout", f"30s timeout, {n_atoms} atoms", elapsed)
    except Exception as e:
        elapsed = time.time() - t0
        return (cid, "exception", f"{type(e).__name__}: {str(e)[:100]}", elapsed)
    finally:
        if tmp_xyz and os.path.exists(tmp_xyz):
            try: os.unlink(tmp_xyz)
            except: pass

def main():
    man = json.load(open("/root/biometaldb-3d/out/full/manifest.json"))
    trex = json.load(open("/root/biometaldb-3d/out/full/trex.json"))
    trex_ids = set(int(k) for k in trex)
    
    tasks = []
    for r in man["records"]:
        cid = r["id"]
        if cid in trex_ids or r.get("status") != "ok":
            continue
        isomers = r.get("isomers", [])
        iso = next((i for i in isomers if i.get("files", {}).get("xyz")), None)
        if not iso:
            continue
        xyz_fn = iso["files"]["xyz"]
        xyz_path = os.path.join("/root/biometaldb-3d/out/full/struct", xyz_fn)
        if not os.path.exists(xyz_path):
            continue
        tasks.append({
            "cid": cid,
            "xyz_path": xyz_path,
            "charge": iso.get("total_charge", 0) or 0,
        })
    
    print(f"Tasks: {len(tasks)}", flush=True)
    
    results = {}
    counts = Counter()
    t0 = time.time()
    
    with ProcessPoolExecutor(max_workers=48) as ex:
        futures = {ex.submit(diag_one, t): t for t in tasks}
        for i, fut in enumerate(as_completed(futures)):
            cid, err_type, detail, elapsed = fut.result()
            results[cid] = (err_type, detail, elapsed)
            counts[err_type] += 1
            if (i + 1) % 100 == 0:
                print(f"  {i+1}/{len(tasks)}: {dict(counts)} ({time.time()-t0:.0f}s)", flush=True)
    
    total_t = time.time() - t0
    print(f"\nDONE ({total_t:.0f}s)")
    print(f"Classification: {dict(counts)}")
    
    # Save detailed results
    with open("/tmp/trex_diag.json", "w") as f:
        json.dump({"counts": dict(counts), "results": {str(k): v for k, v in results.items()}}, f)
    
    # Print samples of each type
    for err_type in sorted(counts.keys()):
        samples = [cid for cid, (t, d, e) in results.items() if t == err_type][:5]
        print(f"\n{err_type} ({counts[err_type]}): samples={samples}")
        if samples:
            cid = samples[0]
            t, d, e = results[cid]
            print(f"  Example #{cid}: {d[:120]} ({e:.1f}s)")

if __name__ == "__main__":
    main()
