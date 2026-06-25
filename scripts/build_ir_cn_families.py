#!/usr/bin/env python3
"""Classify iridium complexes into the two canonical cyclometalated families and
emit a JSON the Flask server renders as a dedicated tab.

Families
--------
  bis_cn_nn  →  [Ir(C^N)2(N^N)]   — two κ²-C,N chelates + one κ²-N,N chelate
                                     (the classic cationic OLED / photoredox motif,
                                      e.g. [Ir(ppy)2(bpy)]+)
  tris_cn    →  [Ir(C^N)3]         — three κ²-C,N chelates
                                     (tris-cyclometalated, e.g. fac-Ir(ppy)3)

Only neutral or cationic complexes are kept (charge_complex >= 0), per request.

Ligand chemistry is read from the precomputed ligand library
(viewer/ligands/ligands.json), which carries per-ligand `denticity` (coordination
number) and `donors` (donor-atom symbols) from the pydentate oracle. A fragment is:
  C^N chelate  ⇔ denticity == 2 and donors == {C, N}
  N^N chelate  ⇔ denticity == 2 and donors == {N, N}
Counterions and solvent (non-coordinating) are ignored; the remaining
coordination sphere must match the target pattern exactly.

Usage:  python scripts/build_ir_cn_families.py \
            [--db data/biometaldb.sqlite] \
            [--ligands viewer/ligands/ligands.json] \
            [--out data/ir_cn_families.json] [--metal Ir]
"""
from __future__ import annotations
import argparse
import json
import os
import sqlite3
from collections import Counter

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NONCOORD_GROUPS = {"Counterion / non-coordinating", "Solvent / aqua"}


def canon(frag: str) -> str:
    """Canonicalize a SMILES fragment the same way ligands_library.py does."""
    m = Chem.MolFromSmiles(frag, sanitize=True) or Chem.MolFromSmiles(frag, sanitize=False)
    return Chem.MolToSmiles(m) if m else frag


def load_ligand_map(path: str) -> dict:
    data = json.load(open(path))
    return {lig["smiles"]: lig for lig in data["ligands"]}


def frag_kind(frag: str, lmap: dict):
    """Return (kind, canonical_smiles) where kind ∈ {CN, NN, noncoord, other, unknown}."""
    c = canon(frag)
    lig = lmap.get(c)
    if lig is None:
        return "unknown", c
    if lig.get("group", "") in NONCOORD_GROUPS:
        return "noncoord", c
    dent = lig.get("denticity")
    donors = tuple(sorted(lig.get("donors") or []))
    if dent == 2 and donors == ("C", "N"):
        return "CN", c
    if dent == 2 and donors == ("N", "N"):
        return "NN", c
    return "other", c


def classify_complex(smiles_ligands: str, lmap: dict):
    """Return ('bis_cn_nn' | 'tris_cn' | None, info_dict)."""
    cn_ligs, nn_ligs, other = [], [], []
    for frag in (f.strip() for f in smiles_ligands.split(".") if f.strip()):
        kind, c = frag_kind(frag, lmap)
        if kind == "CN":
            cn_ligs.append(c)
        elif kind == "NN":
            nn_ligs.append(c)
        elif kind == "noncoord":
            continue  # counterion / solvent — outside the coordination sphere
        else:
            other.append(c)  # unknown / other coordinating ligand → disqualifies

    if other:
        return None, None
    n_cn, n_nn = len(cn_ligs), len(nn_ligs)
    if n_cn == 2 and n_nn == 1:
        return "bis_cn_nn", {"cn_ligands": cn_ligs, "nn_ligand": nn_ligs[0]}
    if n_cn == 3 and n_nn == 0:
        return "tris_cn", {"cn_ligands": cn_ligs}
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=os.path.join(BASE, "data", "biometaldb.sqlite"))
    ap.add_argument("--ligands", default=os.path.join(BASE, "viewer", "ligands", "ligands.json"))
    ap.add_argument("--out", default=os.path.join(BASE, "data", "ir_cn_families.json"))
    ap.add_argument("--metal", default="Ir")
    a = ap.parse_args()

    lmap = load_ligand_map(a.ligands)
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    # measurement counts for display
    meas = dict(con.execute(
        "SELECT complex_id, COUNT(*) FROM measurements GROUP BY complex_id").fetchall())
    rows = con.execute(
        "SELECT id, oxidation_state, charge_complex, smiles_ligands, donor_atoms, has_mol3 "
        "FROM complexes WHERE metal=? ORDER BY id", (a.metal,)).fetchall()
    con.close()

    families = {
        "bis_cn_nn": {
            "key": "bis_cn_nn",
            "formula": "[Ir(C^N)₂(N^N)]",
            "label": "Two C,N-chelates + one N,N-chelate",
            "desc": ("Bis-cyclometalated iridium(III): two κ²-C,N chelating "
                     "ligands and one κ²-N,N chelating ligand "
                     "(e.g. [Ir(ppy)₂(bpy)]⁺). Neutral or cationic."),
            "complexes": [],
        },
        "tris_cn": {
            "key": "tris_cn",
            "formula": "[Ir(C^N)₃]",
            "label": "Three C,N-chelates",
            "desc": ("Tris-cyclometalated iridium(III): three κ²-C,N chelating "
                     "ligands (e.g. fac-Ir(ppy)₃). Neutral or cationic."),
            "complexes": [],
        },
    }

    skipped_anionic = 0
    for cid, ox, charge, smi, donor_atoms, has_mol3 in rows:
        # neutral or cationic only
        if charge is not None and charge < 0:
            fam, _ = classify_complex(smi, lmap)
            if fam:
                skipped_anionic += 1
            continue
        fam, info = classify_complex(smi, lmap)
        if not fam:
            continue
        entry = {
            "id": cid,
            "metal": a.metal,
            "oxidation_state": ox,
            "charge": charge,
            "smiles_ligands": smi,
            "donor_atoms": donor_atoms,
            "has_mol3": bool(has_mol3),
            "n_meas": meas.get(cid, 0),
            **info,
        }
        families[fam]["complexes"].append(entry)

    for f in families.values():
        f["count"] = len(f["complexes"])
        f["charge_breakdown"] = dict(
            Counter(c["charge"] for c in f["complexes"]))

    out = {
        "metal": a.metal,
        "source_db": os.path.basename(a.db),
        "n_bis_cn_nn": families["bis_cn_nn"]["count"],
        "n_tris_cn": families["tris_cn"]["count"],
        "skipped_anionic_matches": skipped_anionic,
        "family_order": ["bis_cn_nn", "tris_cn"],
        "families": families,
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, "w"), ensure_ascii=False)
    print(f"bis_cn_nn (2 C,N + 1 N,N): {out['n_bis_cn_nn']}")
    print(f"tris_cn   (3 C,N):         {out['n_tris_cn']}")
    print(f"skipped anionic matches:   {skipped_anionic}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
