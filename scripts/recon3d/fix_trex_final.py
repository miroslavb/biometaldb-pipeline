"""Final fix: generate ligand-only T-REX for IDs where xyz_to_trex fails.

For complexes where the full xyz_to_trex pipeline crashes (kekulize in
Chem.Mol/RemoveHs, not catchable by SanitizeMol monkey-patch), generate a
simplified T-REX string from SMILES ligands directly:

    Metal{+ox} | L=[ SMILES:..., SMILES:... ]

This preserves the ligand composition information even without MAP.

Usage on hive:
    cd /root/biometaldb-3d
    export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
    ./bin/micromamba run -p ./arch-env python3 polynuclear/fix_trex_final.py
"""
import os, sys, json, time, tempfile, subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter

D_BLOCK = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
           "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
           "Hf","Ta","W","Re","Os","Ir","Pt","Au"}

# Worker: try full xyz_to_trex with deep kekulize patch, fallback to ligand-only
TREX_WORKER_DEEP = """import sys, json, os
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

# Deep monkey-patch: catch kekulize in Mol() constructor AND SanitizeMol AND RemoveHs
_orig_mol_init = Chem.Mol.__init__
_orig_sanitize = Chem.SanitizeMol
_orig_removehs = Chem.RemoveHs

def _patched_mol_init(self, *args, **kwargs):
    try:
        return _orig_mol_init(self, *args, **kwargs)
    except Chem.KekulizeException:
        # Try with sanitize=False
        if 'sanitize' not in kwargs:
            kwargs = dict(kwargs)
            kwargs['sanitize'] = False
            _orig_mol_init(self, *args, **kwargs)
            try:
                from rdkit.Chem import SanitizeFlags
                _orig_sanitize(self, SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE)
            except:
                pass
        else:
            raise

def _patched_sanitize(mol, *args, **kwargs):
    try:
        return _orig_sanitize(mol, *args, **kwargs)
    except Chem.KekulizeException:
        from rdkit.Chem import SanitizeFlags
        flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        return _orig_sanitize(mol, flags)

def _patched_removehs(mol, *args, **kwargs):
    try:
        return _orig_removehs(mol, *args, **kwargs)
    except Chem.KekulizeException:
        # Try sanitizing without kekulize first
        from rdkit.Chem import SanitizeFlags
        try:
            Chem.SanitizeMol(mol, SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE)
        except:
            pass
        return _orig_removehs(mol, *args, **kwargs)

Chem.Mol.__init__ = _patched_mol_init
Chem.SanitizeMol = _patched_sanitize
Chem.RemoveHs = _patched_removehs

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

def build_ligand_only_trex(metal, ox, smiles_ligands, charge, method):
    """Build a simplified T-REX: Metal{+ox} | L=[ SMILES:... ] without MAP."""
    ligands = [s.strip() for s in smiles_ligands.split('.') if s.strip()]
    lig_parts = ", ".join(f"SMILES:{s}" for s in ligands)
    return f"{metal}{{+{ox}}} | L=[ {lig_parts} ]"

def process_one(args_dict):
    cid = args_dict["cid"]
    xyz_path = args_dict["xyz_path"]
    charge = args_dict["charge"]
    method = args_dict["method"]
    metal = args_dict["metal"]
    ox = args_dict["ox"]
    smiles_ligands = args_dict["smiles_ligands"]
    
    n_metals = count_metals(xyz_path)
    if n_metals != 1:
        return (cid, None, "polynuclear")
    
    tmp_xyz = None
    try:
        tmp_xyz = normalize_xyz(xyz_path, charge)
        worker_script = os.path.join(tempfile.gettempdir(), f"_twd_{os.getpid()}.py")
        with open(worker_script, "w") as f:
            f.write(TREX_WORKER_DEEP)
        
        proc = subprocess.run(
            [sys.executable, worker_script, json.dumps({"xyz_path": tmp_xyz, "charge": charge})],
            capture_output=True, text=True, timeout=60,
            env={**os.environ, "PYTHONPATH": "/tmp/trex_repo:" + os.environ.get("PYTHONPATH", "")}
        )
        if proc.returncode == 0 and proc.stdout.strip():
            result = json.loads(proc.stdout.strip())
            return (cid, {"trex": result["trex"], "charge": charge, "method": method}, "ok_full")
        
        # Fallback: ligand-only T-REX
        trex_str = build_ligand_only_trex(metal, ox, smiles_ligands, charge, method)
        return (cid, {"trex": trex_str, "charge": charge, "method": method}, "ok_ligand_only")
    except subprocess.TimeoutExpired:
        # Fallback: ligand-only
        trex_str = build_ligand_only_trex(metal, ox, smiles_ligands, charge, method)
        return (cid, {"trex": trex_str, "charge": charge, "method": method}, "ok_ligand_only")
    except Exception as e:
        # Fallback: ligand-only
        trex_str = build_ligand_only_trex(metal, ox, smiles_ligands, charge, method)
        return (cid, {"trex": trex_str, "charge": charge, "method": method}, "ok_ligand_only")
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
            "method": iso.get("method", ""),
            "metal": r["metal"],
            "ox": r["ox"],
            "smiles_ligands": r.get("smiles_ligands", ""),
        })
    
    print(f"Tasks: {len(tasks)}", flush=True)
    
    n_full = 0
    n_ligand = 0
    n_poly = 0
    t0 = time.time()
    
    with ProcessPoolExecutor(max_workers=48) as ex:
        futures = {ex.submit(process_one, t): t for t in tasks}
        for i, fut in enumerate(as_completed(futures)):
            cid, result, status = fut.result()
            if status == "ok_full":
                trex[str(cid)] = result
                n_full += 1
            elif status == "ok_ligand_only":
                trex[str(cid)] = result
                n_ligand += 1
            elif status == "polynuclear":
                n_poly += 1
            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                print(f"  {i+1}/{len(tasks)}: full={n_full}, ligand_only={n_ligand}, "
                      f"poly={n_poly} ({elapsed:.0f}s)", flush=True)
                with open("/root/biometaldb-3d/out/full/trex.json", "w") as f:
                    json.dump(trex, f, separators=(",", ":"))
    
    with open("/root/biometaldb-3d/out/full/trex.json", "w") as f:
        json.dump(trex, f, separators=(",", ":"))
    
    elapsed = time.time() - t0
    print(f"\nDONE: {n_full} full T-REX, {n_ligand} ligand-only, {n_poly} polynuclear ({elapsed:.0f}s)")
    print(f"Total trex.json: {len(trex)} entries")

if __name__ == "__main__":
    main()
