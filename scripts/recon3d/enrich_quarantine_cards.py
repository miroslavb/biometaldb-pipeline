"""One-time: backfill metadata (smiles, depiction, ligand table, geometry) onto
quarantine/crash records in a manifest that were written as bare stubs by an
older runner. Pulls smiles from the DB, renders the 2D PNG, runs the (crash-safe)
assignment, and overlays it WITHOUT changing the fail verdict. Idempotent.

  python enrich_quarantine_cards.py --db pilot.sqlite --out out/full
"""
from __future__ import annotations
import argparse, json, os, sqlite3
from collections import Counter

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

import assign as A
import coord_oracle as O  # noqa: F401  (oracle cache lives next to the manifest)
from run_ir100 import render_png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", default="out/full")
    args = ap.parse_args()
    man_path = os.path.join(args.out, "manifest.json")
    man = json.load(open(man_path))
    oracle_path = os.path.join(args.out, "oracle.json")
    oracle = json.load(open(oracle_path)) if os.path.exists(oracle_path) else {}

    # bare-stub failures = not built AND no smiles in the record
    stubs = [r for r in man["records"]
             if r.get("status") != "ok" and not r.get("smiles_ligands")]
    if not stubs:
        print("no bare-stub records to enrich")
        return
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    cols = "id, metal, oxidation_state, charge_complex, smiles_ligands, donor_atoms"
    n = 0
    for r in stubs:
        row = con.execute(f"SELECT {cols} FROM complexes WHERE id=?", (r["id"],)).fetchone()
        if not row:
            continue
        cid, metal, ox, charge, smiles, donor_json = row
        r["metal"] = metal
        r["ox"] = ox
        r["charge_complex"] = charge
        r["smiles_ligands"] = smiles
        try:
            png = render_png(cid, metal, ox, smiles,
                             os.path.join(args.out, "img", f"c{cid}.png"))
            if png:
                r["png"] = png
        except Exception:  # noqa: BLE001
            pass
        try:  # assignment is RDKit-only (crash-safe); gives the ligand table
            donor_atoms = json.loads(donor_json) if donor_json else None
            asg = A.assign(metal, ox, smiles, donor_atoms, oracle=oracle)
            r["geometry"] = asg.geometry
            r["cn"] = asg.final_cn
            r["pref_cn"] = asg.pref_cn
            r["assign_conf"] = asg.confidence
            r["hemilabile"] = asg.hemilabile
            r["ligands"] = [{"smiles": u.smiles, "role": u.role, "sites": u.sites,
                             "coordList": u.pocket, "comp": dict(u.pocket_comp),
                             "via": ("oracle" if u.oracle else u.role),
                             "hemilabile": bool(u.oracle and u.oracle.get("hemilabile")),
                             "oracle": u.oracle} for u in asg.units]
        except Exception:  # noqa: BLE001
            pass
        n += 1
    con.close()
    json.dump(man, open(man_path, "w"), indent=2)
    print(f"enriched {n} bare-stub records with smiles/depiction/ligand metadata")


if __name__ == "__main__":
    main()
