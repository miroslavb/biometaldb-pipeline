"""Architector build with a retry ladder + output writers.

Ladder (fast mode): UFF embed -> if no conformer, rescue with GFN2-xTB relax
(handles bulky / low-CN cases like (PPh3)2Au+). xtb mode goes straight to
GFN2-xTB. Returns the lowest-energy conformer as ASE atoms + mol2 connectivity.
"""
from __future__ import annotations
import os
import numpy as np


_RCOV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
         "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39, "As": 1.19, "Se": 1.20,
         "Ru": 1.46, "Rh": 1.42, "Os": 1.44, "Ir": 1.41, "Re": 1.51, "Au": 1.36,
         "Fe": 1.32, "Co": 1.26, "Pt": 1.36, "Pd": 1.39, "Cu": 1.32}


def _stereo_signature(at, metal):
    """Isomer fingerprint = multiset of trans (>150°) donor-element pairs around
    the metal. Distinguishes cis/trans, fac/mer; enantiomers (Λ/Δ) collapse to
    the same signature (mirror images, not separately buildable here)."""
    import numpy as np
    from itertools import combinations
    from collections import Counter
    sym = at.get_chemical_symbols(); xyz = np.asarray(at.get_positions())
    if metal not in sym:
        return None
    mi = sym.index(metal); mr = _RCOV.get(metal, 1.5)
    neigh = [i for i in range(len(sym)) if i != mi
             and np.linalg.norm(xyz[i] - xyz[mi]) <= (mr + _RCOV.get(sym[i], 0.9)) * 1.30]
    trans = []
    for a, b in combinations(neigh, 2):
        v1 = xyz[a] - xyz[mi]; v2 = xyz[b] - xyz[mi]
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
        if np.degrees(np.arccos(np.clip(cos, -1, 1))) > 150:
            trans.append(tuple(sorted((sym[a], sym[b]))))
    return tuple(sorted(Counter(trans).items()))


def _mk_result(conf, method, step):
    import numpy as np
    at = conf["ase_atoms"]
    return {"status": "ok", "method": method, "ladder_step": step,
            "ase_atoms": at, "symbols": at.get_chemical_symbols(),
            "positions": np.asarray(at.get_positions()),
            "mol2string": conf.get("mol2string"), "total_charge": conf.get("total_charge"),
            "n_unpaired": conf.get("calc_n_unpaired_electrons"),
            "energy": conf.get("energy"), "n_atoms": len(at)}


def build_isomers(assignment, mode="xtb", max_isomers=2):
    """Build up to `max_isomers` distinct stereoisomers (cis/trans, fac/mer).

    Enumerates Architector symmetries (n_conformers), deduplicates by stereo
    signature, keeps the lowest-energy representative per distinct isomer.
    Returns a list of result dicts (each with 'isomer' label), or [] on failure.
    """
    from architector import build_complex
    if not assignment.units:
        return []
    if getattr(assignment, "needs_xtb", False):
        mode = "xtb"

    def ligs(with_lt):
        out = []
        for u in assignment.units:
            d = {"smiles": u.smiles, "coordList": list(u.pocket)}
            if with_lt and u.ligType:
                d["ligType"] = u.ligType
            out.append(d)
        return out

    base = dict(metal_ox=assignment.ox, metal_spin=assignment.spin,
                return_only_1=False, n_conformers=4, n_symmetries=12,
                save_init_geos=False, debug=False)
    relax = mode == "xtb"
    amethod = "GFN2-xTB" if relax else "UFF"

    out = None
    for use_lt in (True, False):
        inp = {"core": {"metal": assignment.metal, "coreType": assignment.geometry},
               "ligands": ligs(use_lt),
               "parameters": {**base, "relax": relax, "assemble_method": amethod}}
        try:
            out = build_complex(inp)
        except Exception:
            out = None
        if out:
            break
    if not out:
        return []

    confs = [(k, v) for k, v in out.items()
             if "_init_only" not in k and v.get("ase_atoms") is not None]
    confs.sort(key=lambda kv: kv[1].get("energy", 1e18))
    kept, seen = [], set()
    for _, v in confs:
        sig = _stereo_signature(v["ase_atoms"], assignment.metal)
        if sig in seen:
            continue
        seen.add(sig)
        r = _mk_result(v, amethod + ("+relax" if relax else ""), 0)
        r["stereo_sig"] = str(sig)
        kept.append(r)
        if len(kept) >= max_isomers:
            break
    # label
    if len(kept) == 1:
        kept[0]["isomer"] = "only"
    else:
        for i, r in enumerate(kept):
            r["isomer"] = f"isomer{i+1}"
    return kept


def _best_conf(out):
    if not out:
        return None
    items = [(k, v) for k, v in out.items() if "_init_only" not in k and v.get("ase_atoms") is not None]
    if not items:
        return None
    items.sort(key=lambda kv: kv[1].get("energy", 1e18))
    return items[0][1]


def build(assignment, mode: str = "fast"):
    """mode: 'fast' (UFF, xtb rescue) or 'xtb' (GFN2-xTB relax)."""
    from architector import build_complex

    if not assignment.units:
        return {"status": "no_ligands"}

    # Haptic (Cp/arene) ligands need xtb: UFF leaves the ring floating off-metal.
    if getattr(assignment, "needs_xtb", False) and mode != "xtb":
        mode = "xtb"

    def ligs(with_ligtype):
        out = []
        for u in assignment.units:
            d = {"smiles": u.smiles, "coordList": list(u.pocket)}
            if with_ligtype and u.ligType:
                d["ligType"] = u.ligType
            out.append(d)
        return out

    base = dict(metal_ox=assignment.ox, metal_spin=assignment.spin,
                return_only_1=True, save_init_geos=False, debug=False)

    # (params, use_ligType): ligType first (fast); drop it on retry; xtb rescue.
    if mode == "xtb":
        ladder = [
            (dict(relax=True, assemble_method="GFN2-xTB", n_conformers=2), True),
            (dict(relax=True, assemble_method="GFN2-xTB", n_conformers=2), False),
        ]
    else:
        ladder = [
            (dict(relax=False, assemble_method="UFF", n_conformers=1), True),
            (dict(relax=False, assemble_method="UFF", n_conformers=1), False),
            (dict(relax=True, assemble_method="GFN2-xTB", n_conformers=2), False),  # rescue
        ]

    last_err = None
    for i, (extra, use_lt) in enumerate(ladder):
        inp = {"core": {"metal": assignment.metal, "coreType": assignment.geometry},
               "ligands": ligs(use_lt),
               "parameters": {**base, **extra}}
        try:
            out = build_complex(inp)
        except Exception as e:  # noqa: BLE001
            last_err = f"{type(e).__name__}: {str(e)[:160]}"
            continue
        conf = _best_conf(out)
        if conf is not None:
            at = conf["ase_atoms"]
            return {
                "status": "ok",
                "method": extra["assemble_method"] + ("+relax" if extra["relax"] else ""),
                "ladder_step": i,
                "ase_atoms": at,
                "symbols": at.get_chemical_symbols(),
                "positions": np.asarray(at.get_positions()),
                "mol2string": conf.get("mol2string"),
                "total_charge": conf.get("total_charge"),
                "n_unpaired": conf.get("calc_n_unpaired_electrons"),
                "energy": conf.get("energy"),
                "n_atoms": len(at),
            }
    return {"status": "empty", "error": last_err}


def write_outputs(result, assignment, outdir, cid):
    """Write XYZ (QM input), mol2 (native connectivity incl. metal), SDF."""
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, f"complex_{cid}")
    paths = {}
    at = result["ase_atoms"]
    # XYZ with a charge/mult comment line (consumable by xtb/ORCA wrappers)
    chg = result.get("total_charge")
    mult = (result.get("n_unpaired") or 0) + 1
    xyz = base + ".xyz"
    syms = result["symbols"]; pos = result["positions"]
    with open(xyz, "w") as f:
        f.write(f"{len(syms)}\n")
        f.write(f"cid={cid} metal={assignment.metal} ox={assignment.ox} "
                f"charge={chg} mult={mult} geom={assignment.geometry} "
                f"method={result['method']}\n")
        for s, p in zip(syms, pos):
            f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    paths["xyz"] = xyz
    if result.get("mol2string"):
        m2 = base + ".mol2"
        with open(m2, "w") as f:
            f.write(result["mol2string"])
        paths["mol2"] = m2
        # SDF via OpenBabel (keeps metal-coordination bonds as order-1)
        try:
            from openbabel import openbabel as ob
            conv = ob.OBConversion()
            conv.SetInAndOutFormats("mol2", "sdf")
            obmol = ob.OBMol()
            if conv.ReadString(obmol, result["mol2string"]):
                sdf = base + ".sdf"
                conv.WriteFile(obmol, sdf)
                paths["sdf"] = sdf
        except Exception:
            pass
    return paths
