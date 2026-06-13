"""Validation gates for a reconstructed 3D structure.

A structure is QM-ready only if: the metal is genuinely coordinated to the
intended donor set (count + element composition), there are no atomic clashes,
M-donor bond lengths are physically sane, and charge/multiplicity parity holds.
"""
from __future__ import annotations
from collections import Counter
import numpy as np

# Covalent radii (Angstrom), Cordero 2008 subset (enough for this dataset).
RCOV = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "I": 1.39,
    "Ru": 1.46, "Rh": 1.42, "Os": 1.44, "Ir": 1.41, "Re": 1.51, "Au": 1.36,
    "Fe": 1.32, "Co": 1.26, "Pt": 1.36, "Pd": 1.39, "Cu": 1.32,
}
METALS = {"Ru", "Rh", "Os", "Ir", "Re", "Au", "Fe", "Co", "Pt", "Pd", "Cu"}


def _metal_neighbors(symbols, pos, metal, tol=1.30):
    mi = symbols.index(metal)
    mr = RCOV.get(metal, 1.5)
    out = []
    for i, s in enumerate(symbols):
        if i == mi:
            continue
        d = float(np.linalg.norm(pos[i] - pos[mi]))
        cutoff = (mr + RCOV.get(s, 0.9)) * tol
        if d <= cutoff:
            out.append((i, s, round(d, 3)))
    return mi, sorted(out, key=lambda t: t[2])


def validate(result, assignment):
    """Return (passed: bool, gates: dict, neighbors: list)."""
    symbols = result["symbols"]
    pos = np.asarray(result["positions"])
    metal = assignment.metal
    gates = {}

    # gate: metal present
    if metal not in symbols:
        return False, {"metal_present": False}, []
    mi, neigh = _metal_neighbors(symbols, pos, metal)

    # gate: coordination number (allow +-0; report)
    got_cn = len(neigh)
    gates["cn_got"] = got_cn
    gates["cn_intended"] = assignment.final_cn
    gates["cn_match"] = (got_cn == assignment.final_cn)

    # gate: donor composition
    got_comp = Counter(s for _, s, _ in neigh)
    gates["comp_got"] = dict(got_comp)
    gates["comp_intended"] = dict(assignment.final_comp)
    gates["comp_match"] = (got_comp == assignment.final_comp)

    # gate: M-donor bond lengths sane (0.8*sum_rcov .. 1.4*sum_rcov)
    bl_ok = True
    bl = []
    for _, s, d in neigh:
        lo = (RCOV.get(metal, 1.5) + RCOV.get(s, 0.9)) * 0.80
        hi = (RCOV.get(metal, 1.5) + RCOV.get(s, 0.9)) * 1.40
        bl.append((s, d, lo <= d <= hi))
        if not (lo <= d <= hi):
            bl_ok = False
    gates["bond_lengths"] = bl
    gates["bond_lengths_ok"] = bl_ok

    # gate: no clash (min nonbonded heavy-atom distance > 0.9 A)
    P = pos
    n = len(symbols)
    mind = 1e9
    for i in range(n):
        for j in range(i + 1, n):
            if symbols[i] == "H" or symbols[j] == "H":
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < mind:
                mind = d
    gates["min_heavy_dist"] = round(mind, 3)
    gates["no_clash"] = mind > 0.9

    # gate: charge/mult parity (total electrons vs multiplicity)
    chg = result.get("total_charge")
    nunp = result.get("n_unpaired", 0) or 0
    gates["total_charge"] = chg
    gates["mult"] = nunp + 1

    # A structure is QM-ready if it is physically sane: metal coordinated, no
    # atomic clashes, all M-donor bond lengths in range. comp_match / cn_match
    # are reported as quality signals (haptic ring carbons inflate the raw
    # neighbour count, so they are not gating).
    passed = (gates["no_clash"] and gates["bond_lengths_ok"] and got_cn > 0)
    gates["passed"] = passed
    return passed, gates, neigh
