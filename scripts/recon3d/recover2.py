"""Retry stuck complexes with a conformer-count sweep; emit a mergeable manifest."""
import sys, json, sqlite3, os
sys.path.insert(0, "recon3d")
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
import assign as A, build as B, validate as V, trexio_writer as TX, coord_oracle as O

ids = [int(x) for x in sys.argv[1:]]
OUT = "out/recover2"
os.makedirs(OUT + "/struct", exist_ok=True); os.makedirs(OUT + "/img", exist_ok=True)
con = sqlite3.connect("pilot.sqlite")
def canon(f):
    m = Chem.MolFromSmiles(f, sanitize=True) or Chem.MolFromSmiles(f, sanitize=False)
    return Chem.MolToSmiles(m) if m else f

recs = []
for cid in ids:
    metal, ox, chg, smi, dj = con.execute(
        "SELECT metal,oxidation_state,charge_complex,smiles_ligands,donor_atoms FROM complexes WHERE id=?",
        (cid,)).fetchone()
    da = json.loads(dj) if dj else None
    orc = O.precompute(sorted({canon(f) for f in smi.split(".") if f.strip()}))
    asg = A.assign(metal, ox, smi, da, oracle=orc)
    isos = []
    for nc in (4, 2, 3, 8):
        isos = B.build_isomers(asg, mode="xtb", max_isomers=2, n_conf=nc)
        print(f"#{cid} n_conf={nc} -> {len(isos)} isomers", flush=True)
        if isos:
            break
    rec = {"id": cid, "metal": metal, "ox": ox, "charge_complex": chg, "smiles_ligands": smi,
           "status": "ok" if isos else "no_structure", "geometry": asg.geometry, "cn": asg.final_cn,
           "pref_cn": asg.pref_cn, "spin": asg.spin, "assign_conf": asg.confidence, "flags": asg.flags,
           "hemilabile": asg.hemilabile, "needs_xtb": asg.needs_xtb, "spectators": asg.spectators,
           "ligands": [{"smiles": u.smiles, "role": u.role, "sites": u.sites, "coordList": u.pocket,
                        "comp": dict(u.pocket_comp), "via": ("oracle" if u.oracle else u.role),
                        "hemilabile": bool(u.oracle and u.oracle.get("hemilabile")), "oracle": u.oracle}
                       for u in asg.units], "isomers": []}
    try:
        m = Chem.MolFromSmiles(f"[{metal}+{ox}]." + smi, sanitize=False)
        AllChem.Compute2DCoords(m); Draw.MolToFile(m, f"{OUT}/img/c{cid}.png", size=(520, 380))
        rec["png"] = f"c{cid}.png"
    except Exception:
        rec["png"] = None
    for k, res in enumerate(isos):
        label = res.get("isomer", f"iso{k+1}")
        passed, gates, neigh = V.validate(res, asg)
        paths = B.write_outputs(res, asg, f"{OUT}/struct", f"{cid}_{label}")
        tx = None
        try:
            tx = TX.write_trexio(f"{OUT}/struct/complex_{cid}_{label}.h5", res["symbols"], res["positions"],
                                 res.get("total_charge"), res.get("n_unpaired"), title=f"#{cid} {label}")
        except Exception:
            pass
        rec["isomers"].append({"label": label, "valid": passed, "energy": res.get("energy"),
            "total_charge": res.get("total_charge"), "mult": (res.get("n_unpaired") or 0) + 1,
            "n_atoms": res.get("n_atoms"), "method": res.get("method"), "stereo_sig": res.get("stereo_sig"),
            "metal_neighbors": [[s, d] for _, s, d in neigh],
            "gates": {k2: gates[k2] for k2 in ("no_clash", "bond_lengths_ok", "min_heavy_dist", "cn_got", "comp_got") if k2 in gates},
            "files": {**{k2: os.path.basename(v) for k2, v in paths.items()}, **({"trexio": os.path.basename(tx)} if tx else {})}})
    recs.append(rec)
json.dump({"records": recs}, open(f"{OUT}/manifest.json", "w"), indent=2)
print("DONE", [(r["id"], r["status"], len(r["isomers"])) for r in recs], flush=True)
