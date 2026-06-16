import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, "/root/biometaldb-3d/recon3d")
import sqlite3
import build_poly as BP

con = sqlite3.connect("file:/root/biometaldb-3d/pilot.sqlite?mode=ro", uri=True)
for cid in [9046, 8002, 8868]:
    metal, ox, smi = con.execute(
        "SELECT metal,oxidation_state,smiles_ligands FROM complexes WHERE id=?", (cid,)
    ).fetchone()
    print(f"===== #{cid} {metal}({ox}) =====", flush=True)
    try:
        r = BP.build_polynuclear(metal, ox, smi, uhf=0)
    except Exception as e:
        print(f"  EXC {type(e).__name__}: {str(e)[:160]}", flush=True)
        continue
    if r["status"] != "ok":
        print(f"  FAILED: {r['error']}", flush=True)
        continue
    env = BP.metal_environment(r, metal)
    print(f"  OK seed={r.get('seed')} atoms={r['n_atoms']} q={r['total_charge']} "
          f"seed_min={r.get('seed_min_dist')} E={r['energy']:.1f}", flush=True)
    for k, v in env.items():
        ncp = sum(1 for e in v if e[0] == "C" and 1.9 < e[1] < 2.6)
        print(f"    {k}: n_coord<2.9={len(v)} Cp-C~{ncp} :: {v[:6]}", flush=True)
con.close()
print("ALL_DONE", flush=True)
