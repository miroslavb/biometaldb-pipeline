import warnings; warnings.filterwarnings("ignore")
import sys, os; sys.path.insert(0, "/root/biometaldb-3d/recon3d")
import sqlite3
import build_poly as BP

OUT = "/root/biometaldb-3d/polynuclear/out"
os.makedirs(OUT, exist_ok=True)
IDS = [7323, 7790, 7874, 8002, 8265, 8361, 8382, 8868, 9046]
con = sqlite3.connect("file:/root/biometaldb-3d/pilot.sqlite?mode=ro", uri=True)
ok = 0
for cid in IDS:
    metal, ox, smi = con.execute(
        "SELECT metal,oxidation_state,smiles_ligands FROM complexes WHERE id=?", (cid,)).fetchone()
    print(f"===== #{cid} {metal}({ox}) =====", flush=True)
    try:
        r = BP.build_polynuclear(metal, ox, smi, uhf=0, seeds=(1, 2, 4, 7, 11))
    except Exception as e:
        print(f"  EXC {type(e).__name__}: {str(e)[:140]}", flush=True); continue
    if r.get("status") != "ok":
        print(f"  FAILED: {r.get('error')}", flush=True); continue
    ok += 1
    BP.write_xyz(r, os.path.join(OUT, f"poly_{cid}.xyz"), title=f"BiometalDB #{cid} {metal}({ox})")
    env = BP.metal_environment(r, metal)
    metals_env = " | ".join(f"{k}:{sum(1 for e in v if e[0]=='C' and 1.9<e[1]<2.6)}Cp,"
                            f"{[e for e in v if e[0] in ('S','P','N','O') and e[1]<2.6][:2]}"
                            for k, v in env.items())
    print(f"  OK seed={r.get('seed')} atoms={r['n_atoms']} E={r['energy']:.1f} -> poly_{cid}.xyz | {metals_env}", flush=True)
con.close()
print(f"ALL_DONE built {ok}/{len(IDS)}", flush=True)
