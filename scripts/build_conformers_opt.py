#!/usr/bin/env python3
"""Build zenith_conformers.json from the GFN2-xTB-OPTIMIZED, correct-charge, deduped
conformer ensembles (opt/opt_meta.json), replacing the original UFF-geometry /
charge-0 single-point data. Same output schema the /conformers route consumes.
"""
import json, os, glob

BASE = "/root/.hermes-agent2/biometaldb"
CONFORMER_DIR = "/root/conformer-generation/results/conformers"
COMPOUNDS = "/root/conformer-generation/data/ic50_filtered_compounds.json"
OUT_PATH = os.path.join(BASE, "data", "zenith_conformers.json")
H2KCAL = 627.5094740631

comp = {str(c["id"]): c for c in json.load(open(COMPOUNDS))}
records = []; total_conf = 0; multi = 0; rigid = 0; n_ok = 0; coord_bad = []; failed = []

for d in sorted(glob.glob(os.path.join(CONFORMER_DIR, "complex_*"))):
    cid = int(os.path.basename(d).split("_")[1])
    mpath = os.path.join(d, "opt", "opt_meta.json")
    if not os.path.exists(mpath):
        failed.append(cid); continue
    m = json.load(open(mpath))
    reverted = not m.get("coordination_ok", True)
    orig = {}
    om = os.path.join(d, "metadata.json")
    if os.path.exists(om):
        try: orig = json.load(open(om))
        except Exception: pass
    # conformer_count is authoritative from the SDF actually served by /conformers/sdf
    sdf_path = os.path.join(d, "opt", "conformers_opt_sorted.sdf")
    nk = sum(1 for L in open(sdf_path) if L.strip() == "$$$$") if os.path.exists(sdf_path) else 0
    if nk == 0:
        failed.append(cid); continue
    if reverted:
        # xTB opt broke the coordination sphere -> keep the intact original UFF ensemble
        # and its (charge-0) single-point energies, ordered to match conformers_sorted.sdf.
        e0 = orig.get("gfn2_xtb_energies") or []
        idx = orig.get("gfn2_xtb_sorted_indices") or list(range(len(e0)))
        eh = [e0[i] for i in idx] if e0 else []
        method = "Uniconf UFF (reverted; GFN2-xTB opt broke coordination)"
    else:
        eh = m.get("energies_hartree") or []
        method = "Uniconf seed + GFN2-xTB opt (chrg %s, uhf %s)" % (m.get("chrg"), m.get("uhf"))
    # only attach energies if they line up with the served models, else leave empty (graceful)
    ek = [e * H2KCAL for e in eh] if len(eh) == nk else []
    if not ek:
        eh = []
    nrot = orig.get("n_rotatable_bonds", 0)
    c = comp.get(str(cid), {})
    records.append({
        "id": cid, "conformer_count": nk, "n_rotatable_bonds": nrot, "status": "success",
        "rmsd_max": 0, "reverted": reverted, "method": method,
        "energies_hartree": eh, "energies_kcal": ek,
        "energy_min_kcal": min(ek) if ek else None,
        "energy_max_kcal": max(ek) if ek else None,
        "energy_range_kcal": (max(ek) - min(ek)) if ek else 0,
        "metal": c.get("metal", "Ir"), "oxidation_state": c.get("oxidation_state", 3),
        "charge": c.get("charge", m.get("chrg")), "mw": c.get("mw", 0),
        "smiles_ligands": c.get("smiles_ligands", ""), "min_ic50_dark": c.get("min_ic50_dark"),
        "cn_ligands": c.get("cn_ligands", []), "nn_ligand": c.get("nn_ligand", ""),
        "donor_atoms": c.get("donor_atoms", ""),
        "coordination_ok": m.get("coordination_ok"),
        "has_sdf": os.path.exists(os.path.join(d, "opt", "conformers_opt_sorted.sdf")),
    })
    total_conf += nk; n_ok += 1
    if nk >= 3: multi += 1
    if nk == 1: rigid += 1
    if not m.get("coordination_ok"): coord_bad.append(cid)

records.sort(key=lambda x: x["id"])
out = {
    "total_compounds": n_ok, "total_conformers": total_conf,
    "with_multiple": multi, "rigid_single": rigid,
    "coordination_ok_count": n_ok - len(coord_bad), "coordination_bad": coord_bad,
    "failed_compounds": failed,
    "method_note": "GFN2-xTB optimized (per-complex charge/spin), RMSD+energy deduped distinct minima",
    "records": records,
}
json.dump(out, open(OUT_PATH, "w"), ensure_ascii=False)
print(f"Built {OUT_PATH}: {n_ok} compounds, {total_conf} optimized conformers, "
      f"multi(>=3)={multi}, rigid={rigid}, coord_bad={len(coord_bad)}, failed={len(failed)}")
