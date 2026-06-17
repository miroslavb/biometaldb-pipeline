"""Generate T-REX strings for all built complexes.

For each complex in manifest.json with a valid XYZ structure, generate a T-REX
string using trex.xyz_to_trex. Skip:
  - Polynuclear complexes (multiple metals)
  - Failures

Output: trex.json in the same directory as index.json with
  {id: trex_string} for every successful conversion.
"""
import os
import sys
import json
import time
import shutil
import tempfile
import argparse
import re
from collections import Counter

sys.path.insert(0, '/tmp/trex_repo')
from trex.xyz_to_trex import xyz_to_trex


def normalize_xyz_for_trex(src_path, charge):
    """Create a temporary XYZ file with a T-REX-friendly comment line.

    T-REX's read_xyz_file expects line 1 to look like: 'charge=N ...' and
    errors on '=N' segments. Our generated XYZ files have 'cid=1_only
    metal=Ru ox=2 ...' which contains '=' but no 'charge=' at the right
    position. We rewrite the comment to start with 'charge=N' then append
    a clean label. The file is written to a temp path; caller is responsible
    for cleanup.
    """
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
    D_BLOCK = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
               "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
               "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au"}
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
            if j == 1: continue  # comment
            if not ln.strip(): continue
            sym = ln.split()[0]
            if sym in D_BLOCK: n += 1
            if n_total > 0 and j - 1 >= n_total: break
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/root/biometaldb-3d/out/full")
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--skip-existing", action="store_true")
    args = ap.parse_args()

    manifest_path = args.manifest or os.path.join(args.out, "manifest.json")
    man = json.load(open(manifest_path))

    out_path = os.path.join(args.out, "trex.json")
    trex_data = {}
    if args.skip_existing and os.path.exists(out_path):
        trex_data = json.load(open(out_path))
        print(f"Loaded {len(trex_data)} existing T-REX entries", flush=True)

    recs = man["records"]
    if args.limit:
        recs = recs[:args.limit]

    n_ok = 0
    n_skip_poly = 0
    n_skip_nostruct = 0
    n_fail = 0
    t0 = time.time()
    for i, r in enumerate(recs):
        cid = r["id"]
        if str(cid) in trex_data:
            continue  # already done
        if r.get("status") != "ok":
            n_skip_nostruct += 1
            continue
        # Find first isomer with xyz file
        isomers = r.get("isomers", [])
        isomer = next((i for i in isomers if i.get("files", {}).get("xyz")), None)
        if not isomer:
            n_skip_nostruct += 1
            continue
        xyz_fn = isomer["files"]["xyz"]
        # Try several locations for the XYZ file
        candidates = [
            os.path.join(args.out, "struct", xyz_fn),
            os.path.join(args.out, xyz_fn),
            os.path.join(os.path.dirname(args.out), "out", "struct", xyz_fn),
            os.path.join("/root/biometaldb-3d/out/full/struct", xyz_fn),
        ]
        xyz_path = next((p for p in candidates if os.path.exists(p)), None)
        if not xyz_path:
            n_skip_nostruct += 1
            continue
        # Check polynuclear
        n_metals = count_metals_in_xyz(xyz_path)
        if n_metals != 1:
            n_skip_poly += 1
            continue
        # Generate T-REX (via normalized xyz with charge= comment)
        chg = isomer.get("total_charge", 0) or 0
        tmp_xyz = None
        try:
            tmp_xyz = normalize_xyz_for_trex(xyz_path, chg)
            trex_str = xyz_to_trex(tmp_xyz, overall_charge=chg)
            trex_data[str(cid)] = {
                "trex": trex_str,
                "charge": chg,
                "method": isomer.get("method", ""),
            }
            n_ok += 1
        except Exception as e:
            n_fail += 1
            if n_fail <= 3:
                print(f"  FAIL #{cid}: {type(e).__name__}: {str(e)[:80]}", flush=True)
        finally:
            if tmp_xyz and os.path.exists(tmp_xyz):
                try: os.unlink(tmp_xyz)
                except: pass
        if (i + 1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"  {i+1}/{len(recs)}: ok={n_ok}, poly_skip={n_skip_poly}, "
                  f"no_struct={n_skip_nostruct}, fail={n_fail} ({elapsed:.0f}s)",
                  flush=True)
            # save progress
            with open(out_path, "w") as f:
                json.dump(trex_data, f, separators=(",", ":"))

    with open(out_path, "w") as f:
        json.dump(trex_data, f, separators=(",", ":"))
    elapsed = time.time() - t0
    print(f"DONE: {n_ok} T-REX generated, {n_skip_poly} polynuclear skipped, "
          f"{n_fail} failed, {n_skip_nostruct} no_struct ({elapsed:.0f}s)")


if __name__ == "__main__":
    main()
