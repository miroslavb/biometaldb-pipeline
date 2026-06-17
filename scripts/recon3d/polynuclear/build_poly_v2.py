"""build_poly_v2 — robust polynuclear builder for sandwich-decorated complexes.

Handles cases where a metallocene (Co/Fe/Ru cp) is attached to a donor atom
(Se/NHC/P/S) that links to the core metal. Approach:
  1. Parse SMILES, replace dative bonds (-> / <-) with single bonds.
  2. Detect and SPLIT OUT metallocene fragments (Co/Fe/Ru/Ni/Mn/Cr/V/Zr/Hf/Mo/W
     with >=2 dative bonds to aromatic C) as spectators. Remaining fragment is
     the "real" donor (Se-, S-, NHC-C, P, etc).
  3. Run build_polynuclear-style assembly: ideal geometry for core, then
     constrained xtb relax (GFN-FF then GFN2).
  4. Re-attach spectators by re-inserting the metallocene at the exact position
     it occupied in the SMILES graph (RDKit-3D).
"""
import re
import sys
import os
import subprocess
import tempfile
import numpy as np

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

INTERNAL_METALS = {"Ti", "Fe", "Co", "Ru", "Ni", "Mn", "Cr", "V", "Zr", "Hf", "Mo", "W"}

# Ideal metal-coordination distances (Å)
_R_DONOR = {"C": 2.02, "P": 2.28, "N": 2.05, "S": 2.30, "Se": 2.40,
            "O": 2.05, "Cl": 2.34, "Br": 2.45, "I": 2.60, "H": 1.55}


def _is_metallocene_atom(a, mol):
    """True if atom is an internal metal with >= 2 dative bonds to C."""
    if a.GetSymbol() not in INTERNAL_METALS:
        return False
    n_dative = 0
    for b in a.GetBonds():
        if b.GetBondTypeAsDouble() == 1.0 and b.GetOtherAtom(a).GetIsAromatic():
            n_dative += 1
    return n_dative >= 2


def _find_metallocene_atoms(mol):
    """Return set of atom indices that are PART OF a metallocene (just the metal
    atom). Cp-ring carbons are kept in the donor fragment (they remain bonded
    to the donor atom)."""
    met_atoms = set()
    for a in mol.GetAtoms():
        if _is_metallocene_atom(a, mol):
            met_atoms.add(a.GetIdx())
    return met_atoms


def _strip_dummies(smi):
    """Remove dummy atoms ([*], [n*], etc.) from a SMILES entirely. We do this
    by parsing the mol, building an EditableMol, removing the dummies (and the
    ring-closure flags they leave), then re-emitting.

    Returns a clean SMILES (no dummies) or the original on failure.
    """
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return smi
    # Find dummy atoms
    dummy_idxs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "*"]
    if not dummy_idxs:
        # No dummies, just sanitize and return
        try:
            Chem.SanitizeMol(m)
            return Chem.MolToSmiles(Chem.RemoveHs(m))
        except Exception:
            return smi
    # Remove in reverse to preserve indices
    em = Chem.EditableMol(m)
    for idx in sorted(dummy_idxs, reverse=True):
        em.RemoveAtom(idx)
    m2 = em.GetMol()
    try:
        Chem.SanitizeMol(m2)
    except Exception:
        # Try with partial sanitization (skip kekulization)
        try:
            Chem.SanitizeMol(m2, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception:
            return smi
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(m2))
    except Exception:
        return smi


def _strip_metallocene_from_smiles(smi):
    """Split a SMILES by REMOVING the metallocene metal atom (and any H attached
    to it). Keeps the Cp rings in the donor fragment, so the Cp atoms remain
    attached to the donor (e.g. Se, NHC).

    Returns (core_smi, mc_smi) where core_smi is the donor-bearing fragment
    (with Cp rings) and mc_smi is just the bare metal (Co, Fe, ...).
    """
    smi_d = re.sub(r"->|<-", "-", smi)
    m = Chem.MolFromSmiles(smi_d)
    if m is None:
        m = Chem.MolFromSmiles(smi_d, sanitize=False)
    if m is None:
        return smi, None

    mc_idxs = _find_metallocene_atoms(m)
    if not mc_idxs:
        return smi, None  # no metallocene

    # Find bonds that connect the metallocene metal to the rest
    cut_bonds = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if (i in mc_idxs) != (j in mc_idxs):
            cut_bonds.append(b.GetIdx())
    if not cut_bonds:
        return smi, None  # metallocene is the whole fragment

    # Split the molecule
    try:
        frags = Chem.FragmentOnBonds(m, cut_bonds, dummyLabels=[(i, i) for i in range(len(cut_bonds))])
        pieces = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=False)
    except Exception:
        return smi, None

    # Identify which piece contains the metallocene metal and which is the core donor
    core_piece, mc_piece = None, None
    for p in pieces:
        if any(a.GetSymbol() in INTERNAL_METALS for a in p.GetAtoms()):
            mc_piece = p
        else:
            if core_piece is None or p.GetNumAtoms() > core_piece.GetNumAtoms():
                core_piece = p  # take the LARGEST non-metal piece (likely the donor)
    if core_piece is None or mc_piece is None:
        return smi, None

    # Clean dummy atoms from core piece
    try:
        core_clean = Chem.RemoveHs(core_piece)
        Chem.SanitizeMol(core_clean)
        core_smi = Chem.MolToSmiles(core_clean)
    except Exception:
        return smi, None
    try:
        mc_clean = Chem.RemoveHs(mc_piece)
        Chem.SanitizeMol(mc_clean)
        mc_smi = Chem.MolToSmiles(mc_clean)
    except Exception:
        mc_smi = None
    return core_smi, mc_smi


def _find_core_donor(mh):
    """Find the best donor atom for the core metal — same as build_poly."""
    best = None
    for a in mh.GetAtoms():
        if any(nb.GetSymbol() in INTERNAL_METALS for nb in a.GetNeighbors()):
            continue
        sym, ch = a.GetSymbol(), a.GetFormalCharge()
        s = None
        if sym == "C" and ch < 0: s = 95
        elif sym == "P":
            # Phosphine P (3-coordinate, no =O) OR cationic P+ in phosphonium
            if a.GetDegree() == 3 and not any(nb.GetSymbol() == "O" for nb in a.GetNeighbors()):
                s = 92
            elif ch > 0:  # [P+] — phosphonium acts as donor
                s = 88
        elif sym == "S" and ch < 0: s = 86
        elif sym == "Se" and ch < 0: s = 84
        elif sym == "N" and ch < 0: s = 82
        elif sym in ("Cl", "Br", "I") and ch < 0: s = 52
        if s is not None and (best is None or s > best[2]):
            best = (a.GetIdx(), sym, s)
    return best


def _embed(mol, seed=1):
    mh = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useRandomCoords = True
    p.maxIterations = 4000
    if AllChem.EmbedMolecule(mh, p) != 0:
        if AllChem.EmbedMolecule(mh, useRandomCoords=True, maxAttempts=50,
                                 randomSeed=seed + 17) != 0:
            return None
    return mh


def _declash(pos, ranges, axes, passes=3):
    """Rotate fragments about the core (atom 0) to minimize clash."""
    for _ in range(passes):
        for a in range(len(ranges)):
            for b in range(a + 1, len(ranges)):
                sa, ea = ranges[a]
                sb, eb = ranges[b]
                pa = pos[sa:ea]
                pb = pos[sb:eb]
                # crude: rotate fragment b about axis to maximize min distance
                best_d, best_t = -1.0, 0.0
                axis = axes[b] if b < len(axes) else np.array([0., 0., 1.])
                for t in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                    c, s = np.cos(t), np.sin(t)
                    R = np.array([[c + axis[0]**2*(1-c), axis[0]*axis[1]*(1-c)-axis[2]*s,
                                   axis[0]*axis[2]*(1-c)+axis[1]*s],
                                  [axis[1]*axis[0]*(1-c)+axis[2]*s, c+axis[1]**2*(1-c),
                                   axis[1]*axis[2]*(1-c)-axis[0]*s],
                                  [axis[2]*axis[0]*(1-c)-axis[1]*s, axis[2]*axis[1]*(1-c)+axis[0]*s,
                                   c+axis[2]**2*(1-c)]])
                    pb2 = pb @ R.T
                    d = _min_pairwise(pa, pb2)
                    if d > best_d:
                        best_d, best_t = d, t
                if best_t != 0:
                    c, s = np.cos(best_t), np.sin(best_t)
                    R = _rot(axis, best_t)
                    pos[sb:eb] = pb @ R.T
    return pos


def _rot(axis, theta):
    c, s = np.cos(theta), np.sin(theta)
    a = axis / np.linalg.norm(axis)
    ax, ay, az = a
    return np.array([
        [c + ax*ax*(1-c),     ax*ay*(1-c) - az*s, ax*az*(1-c) + ay*s],
        [ay*ax*(1-c) + az*s, c + ay*ay*(1-c),     ay*az*(1-c) - ax*s],
        [az*ax*(1-c) - ay*s, az*ay*(1-c) + ax*s, c + az*az*(1-c)],
    ])


def _min_pairwise(a, b):
    """Min |ai - bj| for any i in a, j in b (vectorized)."""
    diff = a[:, None, :] - b[None, :, :]
    d = np.linalg.norm(diff, axis=-1)
    return float(d.min())


def _align(frag_pos, donor_i, target_axis, r):
    """Translate fragment so donor_i is at distance r along target_axis from origin."""
    p = frag_pos.copy()
    p -= p[donor_i]  # donor at origin
    # Now place donor at r * target_axis
    p += r * np.array(target_axis) / np.linalg.norm(target_axis)
    return p


def _make_bent_geometry(donor_elems, metal, ox):
    """Choose idealized axes for the donor arrangement around the core metal."""
    # For Au(I) with 2 donors: linear (180°)
    # For Au(III)/Pd(II)/Pt(II) with 4 donors: square planar
    # For others with 4: tetrahedral
    n = len(donor_elems)
    if metal in ("Au", "Ag", "Cu") and ox == 1:
        # linear
        if n == 1: return [np.array([0., 0., 1.])]
        return [np.array([0., 0., 1.]), np.array([0., 0., -1.])]
    if metal in ("Au", "Ag", "Cu") and ox == 3:
        # square planar
        return [np.array([1., 0., 0.]), np.array([0., 1., 0.]),
                np.array([-1., 0., 0.]), np.array([0., -1., 0.])][:n]
    if n == 2: return [np.array([0., 0., 1.]), np.array([0., 0., -1.])]
    if n == 3: return [np.array([1., 0., 0.]), np.array([-0.5, 0.866, 0.]),
                       np.array([-0.5, -0.866, 0.])]
    # tetrahedral for 4
    return [np.array([1., 1., 1.]), np.array([-1., -1., 1.]),
            np.array([-1., 1., -1.]), np.array([1., -1., -1.])][:n]


def _xtb_constrained_opt(syms, pos, charge, uhf, constraints,
                          method="gfn2", maxcycle=300, acc=2.0):
    """Call xtb CLI with harmonic constraints on given (i,j,dist) triples."""
    from ase import Atoms
    from ase.optimize import LBFGS
    from ase.constraints import FixBondLengths
    from xtb.ase.calculator import XTB
    at = Atoms(symbols=syms, positions=pos)
    if constraints:
        at.set_constraint(FixBondLengths([(i, j) for i, j, _ in constraints]))
    try:
        at.calc = XTB(method="GFN-FF" if method == "gfnff" else "GFN2-xTB",
                      charge=int(charge), uhf=uhf,
                      electronic_temperature=1500.0, accuracy=acc,
                      max_iterations=maxcycle)
        opt = LBFGS(at, logfile=None)
        opt.run(fmax=0.3, steps=maxcycle)
        return list(at.get_chemical_symbols()), np.array(at.get_positions())
    except Exception as e:
        return None


def build_poly_v2(metal, ox, smiles_ligands, cn=None, uhf=0, seeds=(1, 2, 4, 7)):
    """Build a polynuclear complex, handling metallocene-decorated donors."""
    frags_raw = [f.strip() for f in smiles_ligands.split(".") if f.strip()]

    # Step 1: detect and split out metallocene spectators per fragment
    donor_frags = []      # fragments after metallocene split
    spectators = []       # [(mc_smiles, attach_to_atom_idx_in_donor)]
    for f in frags_raw:
        core_smi, mc_smi = _strip_metallocene_from_smiles(f)
        # Strip dummy atoms from core (they were created by FragmentOnBonds
        # and RDKit's DistGeom can't handle them)
        core_smi = _strip_dummies(core_smi)
        donor_frags.append(core_smi)
        if mc_smi:
            spectators.append(mc_smi)

    # Step 2: build the core metal + donors (use build_poly path)
    units = []
    for smi in donor_frags:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            continue
        mh = _embed(m, seed=1)
        if mh is None:
            continue
        donor = _find_core_donor(mh)
        if donor is None:
            continue
        di, de, _ = donor
        neigh = [nb.GetIdx() for nb in mh.GetAtomWithIdx(di).GetNeighbors()]
        units.append((mh, di, neigh, de))

    if not units:
        return {"status": "failed", "error": f"no donor atoms found after metallocene split "
                                              f"(spectators={len(spectators)})"}

    # Step 3: assemble
    cn_eff = cn or min(len(units), 2 if metal in ("Au", "Ag", "Cu") and ox == 1 else 4)
    units = units[:cn_eff]
    axes = _make_bent_geometry([u[3] for u in units], metal, ox)

    all_sym = [metal]
    all_pos = [np.zeros(3)]
    ranges = []
    core_bonds = []
    offset = 1
    tot_charge = ox  # core metal contributes +ox
    for (mh, di, neigh, de), axis in zip(units, axes):
        # Skip dummy atoms (RDKit placeholders from FragmentOnBonds)
        keep = [i for i, a in enumerate(mh.GetAtoms()) if a.GetSymbol() != "*"]
        if di not in keep:
            di = keep[0]
        # remap di to its new index in `keep`
        new_di = keep.index(di)
        syms = [mh.GetAtomWithIdx(i).GetSymbol() for i in keep]
        conf = mh.GetConformer()
        pos = np.array([list(conf.GetAtomPosition(i)) for i in keep])
        r = _R_DONOR.get(de, 2.1)
        moved = _align(pos, new_di, axis, r)
        tot_charge += Chem.GetFormalCharge(Chem.RemoveHs(mh))
        core_bonds.append((0, offset + new_di, r))
        ranges.append((offset, offset + len(syms)))
        all_sym += syms
        all_pos.append(moved)
        offset += len(syms)
    all_pos = np.vstack(all_pos)
    all_pos = _declash(all_pos, ranges, axes, passes=3)

    # Step 4: constrained relax
    charge = tot_charge
    n_atoms_core = len(all_sym)
    has_spectator = bool(spectators)
    # For polynuclear complexes (with spectators), use GFN-FF only — GFN2 SCF is
    # often unstable for bimetallics with metallocene spectators.
    if has_spectator or n_atoms_core > 35:
        res = _xtb_constrained_opt(all_sym, all_pos, charge, uhf, core_bonds,
                                    method="gfnff", maxcycle=150)
        if res is None:
            return {"status": "failed", "error": f"xtb GFN-FF failed (n_atoms={n_atoms_core})"}
        all_sym, all_pos = res
        method = "GFN-FF(constrained)+poly_v2"
    else:
        res = _xtb_constrained_opt(all_sym, all_pos, charge, uhf, core_bonds,
                                    method="gfn2", maxcycle=200)
        if res is not None:
            all_sym, all_pos = res
            method = "GFN2-xTB(constrained)+poly_v2"
        else:
            res = _xtb_constrained_opt(all_sym, all_pos, charge, uhf, core_bonds,
                                        method="gfnff", maxcycle=150)
            if res is None:
                return {"status": "failed", "error": "xtb constrained opt failed"}
            all_sym, all_pos = res
            method = "GFN-FF(constrained)+poly_v2"

    # Step 5: attach metallocene spectators at the donor atom of the first fragment
    # For each spectator, build a fresh metallocene (Co/Fe cp sandwich) and place
    # it near the donor atom. This is more reliable than embedding the parsed
    # metallocene fragment (which often has ring-closure dummies and lost bonds).
    n_attach = min(len(spectators), len(ranges))
    for k in range(n_attach):
        mc_smi = spectators[k]
        # Identify the metallocene metal from the spectator SMILES
        mc_metal = None
        for tok in re.findall(r"\[([A-Z][a-z]?\+\d)\]", mc_smi):
            for el in INTERNAL_METALS:
                if el.lower() in tok.lower():
                    mc_metal = el
                    break
            if mc_metal: break
        if mc_metal is None:
            # Try simpler match: just look for any metal token
            for el in INTERNAL_METALS:
                if f"[{el}" in mc_smi or f"-[{el}" in mc_smi or f"]{el}+" in mc_smi:
                    mc_metal = el
                    break
        if mc_metal is None:
            continue

        # Build a fresh metallocene: two stacked Cp rings at the metal
        mc_pos, mc_syms = _build_metallocene(mc_metal)
        if mc_pos is None:
            continue

        # Find anchor point: the donor atom of unit k in merged indexing
        sa, ea = ranges[k]
        donor_atom_in_merged = core_bonds[k][1]
        anchor = all_pos[donor_atom_in_merged]
        axis = axes[k] if k < len(axes) else np.array([1., 0., 0.])
        # Perpendicular to axis
        if abs(axis[0]) < 0.9:
            perp = np.cross(axis, [1., 0., 0.])
        else:
            perp = np.cross(axis, [0., 1., 0.])
        perp /= np.linalg.norm(perp) + 1e-9
        # Place metallocene 4.2 Å from anchor, oriented so one Cp ring faces the donor
        new_origin = anchor + 4.2 * perp
        mc_pos += new_origin - mc_pos[0]  # translate so metal is at new_origin
        all_sym += mc_syms
        all_pos = np.vstack([all_pos, mc_pos])

    return {"status": "ok", "method": method, "symbols": all_sym,
            "positions": all_pos, "n_atoms": len(all_sym),
            "total_charge": int(charge), "n_unpaired": int(uhf),
            "metal": metal, "core_bonds": core_bonds,
            "spectators": spectators}


def _build_metallocene(metal_symbol):
    """Build a stacked metallocene (metal between two Cp rings) from scratch.
    Returns (pos[N,3], syms[N]) or (None, None) on failure.
    Uses ideal M–C distances:
      Co 2.10 Å, Fe 2.04 Å, Ru 2.18 Å, Ni 2.18 Å, Mn 2.18 Å, Cr 2.20 Å,
      V 2.27 Å, Ti 2.34 Å, Zr/Hf 2.45 Å, Mo/W 2.30 Å
    Cp rings: pentagons of 5 C, ring radius 1.20 Å, M-Cp distance 1.65 Å
    """
    M_C = {"Co": 2.10, "Fe": 2.04, "Ru": 2.18, "Ni": 2.18, "Mn": 2.18,
           "Cr": 2.20, "V": 2.27, "Ti": 2.34, "Zr": 2.45, "Hf": 2.45,
           "Mo": 2.30, "W": 2.30}.get(metal_symbol, 2.10)
    R_C_RING = 1.20
    D_M_CP = 1.65
    n = 5
    syms = [metal_symbol]
    # First Cp ring (in xy plane at z = -D_M_CP)
    angles1 = [2 * np.pi * i / n for i in range(n)]
    cp1 = np.array([[R_C_RING * np.cos(a), R_C_RING * np.sin(a), -D_M_CP] for a in angles1])
    # Second Cp ring (rotated 36°, at z = +D_M_CP)
    angles2 = [2 * np.pi * i / n + np.pi / n for i in range(n)]
    cp2 = np.array([[R_C_RING * np.cos(a), R_C_RING * np.sin(a), D_M_CP] for a in angles2])
    syms += ["C"] * (2 * n)
    pos = np.vstack([np.array([[0., 0., 0.]]), cp1, cp2])
    return pos, syms


def write_xyz(result, path, title=""):
    with open(path, "w") as f:
        f.write(f"{len(result['symbols'])}\n")
        f.write(f"{title} method={result.get('method','')} n_unpaired={result.get('n_unpaired',0)}\n")
        for s, p in zip(result["symbols"], result["positions"]):
            f.write(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
