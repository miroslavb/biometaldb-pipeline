"""Fix T-REX failures: retry with kekulize bypass + extended timeout.

Strategy:
  - KekulizeException IDs: monkey-patch Chem.SanitizeMol to skip SANITIZE_KEKULIZE
  - Timeout IDs: extend timeout to 60s
  - ok_with_30s IDs: just re-run (they work with longer timeout)
  - exit_with_warnings: try with kekulize bypass

Usage on hive:
    cd /root/biometaldb-3d
    export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
    ./bin/micromamba run -p ./arch-env python3 polynuclear/fix_trex_fails.py
"""
import os, sys, json, time, tempfile, subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter

D_BLOCK = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
           "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
           "Hf","Ta","W","Re","Os","Ir","Pt","Au"}

# Worker script that monkey-patches SanitizeMol to skip kekulize
TREX_WORKER_PATCHED = """import sys, json, os
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

# Monkey-patch: skip kekulize during sanitization
_orig_sanitize = Chem.SanitizeMol
def _patched_sanitize(mol, *args, **kwargs):
    try:
        return _orig_sanitize(mol, *args, **kwargs)
    except Chem.KekulizeException:
        # Retry without kekulize
        from rdkit.Chem import SanitizeFlags
        flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        return _orig_sanitize(mol, flags)
Chem.SanitizeMol = _patched_sanitize

sys.path.insert(0, '/tmp/trex_repo')
from trex.xyz_to_trex import xyz_to_trex
args = json.loads(sys.argv[1])
result = xyz_to_trex(args['xyz_path'], overall_charge=args['charge'])
print(json.dumps({'ok': True, 'trex': result}))
"""

TREX_WORKER_STANDARD = """import sys, json
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
            if n_total is not None and n_total > 0 and j - 1 >= n_total: break
    return n

def process_one(args_dict):
    cid = args_dict["cid"]
    xyz_path = args_dict["xyz_path"]
    charge = args_dict["charge"]
    method = args_dict["method"]
    use_patched = args_dict.get("patched", False)
    timeout = args_dict.get("timeout", 30)
    
    n_metals = count_metals(xyz_path)
    if n_metals != 1:
        return (cid, None, "polynuclear")
    
    tmp_xyz = None
    try:
        tmp_xyz = normalize_xyz(xyz_path, charge)
        worker_code = TREX_WORKER_PATCHED if use_patched else TREX_WORKER_STANDARD
        worker_script = os.path.join(tempfile.gettempdir(), f"_twf_{os.getpid()}.py")
        with open(worker_script, "w") as f:
            f.write(worker_code)
        
        proc = subprocess.run(
            [sys.executable, worker_script, json.dumps({"xyz_path": tmp_xyz, "charge": charge})],
            capture_output=True, text=True, timeout=timeout,
            env={**os.environ, "PYTHONPATH": "/tmp/trex_repo:" + os.environ.get("PYTHONPATH", "")}
        )
        if proc.returncode != 0:
            return (cid, None, f"exit={proc.returncode}: {proc.stderr[:120]}")
        result = json.loads(proc.stdout.strip())
        return (cid, {"trex": result["trex"], "charge": charge, "method": method}, "ok")
    except subprocess.TimeoutExpired:
        return (cid, None, f"timeout_{timeout}s")
    except Exception as e:
        return (cid, None, f"{type(e).__name__}: {str(e)[:100]}")
    finally:
        if tmp_xyz and os.path.exists(tmp_xyz):
            try: os.unlink(tmp_xyz)
            except: pass

def main():
    man = json.load(open("/root/biometaldb-3d/out/full/manifest.json"))
    trex = json.load(open("/root/biometaldb-3d/out/full/trex.json"))
    trex_ids = set(int(k) for k in trex)
    
    # Load diagnosis results
    diag = json.load(open("/tmp/trex_diag.json"))
    diag_results = diag["results"]
    
    tasks = []
    stats = Counter()
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
        
        err_type = diag_results.get(str(cid), ("unknown", "", 0))[0]
        
        if err_type == "polynuclear":
            stats["skip_poly"] += 1
            continue
        
        task = {
            "cid": cid,
            "xyz_path": xyz_path,
            "charge": iso.get("total_charge", 0) or 0,
            "method": iso.get("method", ""),
        }
        
        if err_type in ("ok",):
            # Already works with 30s, just re-run
            task["timeout"] = 30
            task["patched"] = False
            stats["retry_ok"] += 1
        elif err_type in ("other_exit", "exception", "exit_with_warnings"):
            # Kekulize-related — use patched sanitize
            task["timeout"] = 30
            task["patched"] = True
            stats["retry_patched"] += 1
        elif err_type == "timeout":
            # Slow — extend timeout
            task["timeout"] = 60
            task["patched"] = False
            stats["retry_timeout60"] += 1
        else:
            task["timeout"] = 30
            task["patched"] = True
            stats["retry_unknown"] += 1
        
        tasks.append(task)
    
    print(f"Tasks: {len(tasks)} — {dict(stats)}", flush=True)
    
    n_ok = 0
    n_fail = 0
    fail_types = Counter()
    t0 = time.time()
    
    with ProcessPoolExecutor(max_workers=48) as ex:
        futures = {ex.submit(process_one, t): t for t in tasks}
        for i, fut in enumerate(as_completed(futures)):
            cid, result, status = fut.result()
            if status == "ok":
                trex[str(cid)] = result
                n_ok += 1
            else:
                n_fail += 1
                fail_types[status] += 1
            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                print(f"  {i+1}/{len(tasks)}: ok={n_ok}, fail={n_fail} ({elapsed:.0f}s)", flush=True)
                # Save progress
                with open("/root/biometaldb-3d/out/full/trex.json", "w") as f:
                    json.dump(trex, f, separators=(",", ":"))
    
    # Final save
    with open("/root/biometaldb-3d/out/full/trex.json", "w") as f:
        json.dump(trex, f, separators=(",", ":"))
    
    elapsed = time.time() - t0
    print(f"\nDONE: {n_ok} recovered, {n_fail} still failed ({elapsed:.0f}s)")
    print(f"Fail types: {dict(fail_types)}")
    print(f"Total trex.json: {len(trex)} entries")

if __name__ == "__main__":
    main()
