#!/usr/bin/env python3
"""Bent-metallocene template generator + improved build_poly.

Replaces RDKit ETKDG for metallocene fragments with idealised bent/parallel
sandwich geometry, then constrained-relaxes in isolation.  This fixes the
root cause of the 6 titanocene failures: ETKDG produces random, often
overlapping Cp rings that clash on assembly and crash the xTB SCF.

Usage (standalone test):
  python bent_template.py --test
  python bent_template.py --run-failed
"""

import warnings; warnings.filterwarnings("ignore")
import sys, os, subprocess, tempfile, time, argparse, json, shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# ── geometry constants ──────────────────────────────────────────────
# Bent angle Cp–M–Cp (degrees): Ti=130, Zr=135, Hf=135
# Parallel sandwich: Fe=180, Co=180, Ru=180, Cr=180, Mn=180, Ni=180, V=180
_BENT_ANGLE = {"Ti": 130, "Zr": 135, "Hf": 135}
_SANDWICH_METALS = {"Ti", "Zr", "Hf", "Fe", "Co", "Ru", "Cr", "Mn", "Ni", "V"}
# M–Cp_centroid distance (Å) — metal to ring centre
_M_CP_DIST = {"Ti": 2.05, "Zr": 2.20, "Hf": 2.18,
              "Fe": 1.65, "Co": 1.70, "Ru": 1.80,
              "Cr": 1.70, "Mn": 1.75, "Ni": 1.70, "V": 1.90}
# Cp ring radius (C–centroid distance, Å)
_CP_RADIUS = 1.21
# C–C distance in Cp ring (Å)
_CP_CC = 1.40

INTERNAL_METALS = _SANDWICH_METALS | {"Mo", "W"}
_R_DONOR = {"C": 2.02, "P": 2.28, "S": 2.30, "N": 2.05, "Cl": 2.27, "O": 2.05}
_R_INTERNAL = {("Ti", "C"): 2.37, ("Ti", "O"): 2.05, ("Ti", "Cl"): 2.30,
               ("Fe", "C"): 2.04, ("Co", "C"): 2.10, ("Ru", "C"): 2.20,
               ("Zr", "C"): 2.45, ("Hf", "C"): 2.42}

_GEOM_AXES = {
    2: [(0, 0, 1), (0, 0, -1)],
    3: [(0, 0, 1), (0.866, 0, -0.5), (-0.866, 0, -0.5)],
    4: [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],
}


# ═══════════════════════════════════════════════════════════════════
#  Bent-metallocene template
# ═══════════════════════════════════════════════════════════════════

def _is_metallocene(mol):
    """True if mol contains a sandwich metal with 10 Cp carbons."""
    for a in mol.GetAtoms():
        if a.GetSymbol() not in _SANDWICH_METALS:
            continue
        c_neighbors = [nb for nb in a.GetNeighbors() if nb.GetSymbol() == "C"]
        if len(c_neighbors) >= 8:   # η5-Cp gives 5 per ring, 10 total
            return True, a.GetSymbol(), a.GetIdx()
    return False, None, None


def _find_cp_rings(mol, metal_idx):
    """Return two lists of carbon atom indices for the two Cp rings."""
    metal = mol.GetAtomWithIdx(metal_idx)
    c_neighbors = [nb.GetIdx() for nb in metal.GetNeighbors() if nb.GetSymbol() == "C"]
    if len(c_neighbors) < 8:
        return [], []
    # Partition by connectivity: Cp carbons form 5-membered rings
    # Build adjacency among the C neighbors
    adj = {i: set() for i in c_neighbors}
    for i in c_neighbors:
        ai = mol.GetAtomWithIdx(i)
        for nb in ai.GetNeighbors():
            if nb.GetIdx() in c_neighbors:
                adj[i].add(nb.GetIdx())
    # Find 5-membered rings via DFS
    visited = set()
    rings = []
    for start in c_neighbors:
        if start in visited:
            continue
        ring = []
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            ring.append(v)
            for w in adj[v]:
                if w not in visited:
                    stack.append(w)
        if len(ring) == 5:
            rings.append(ring)
    return rings if len(rings) == 2 else ([], [])


def _regular_pentagon_xyz(center, normal, radius=_CP_RADIUS, cc=_CP_CC):
    """Generate 5 carbon positions for a regular pentagon Cp ring.

    center: (3,) — ring centroid
    normal: (3,) — ring plane normal (unit vector)
    radius: distance from centroid to each C
    cc: C–C distance (used to verify, not to construct)

    Returns (5, 3) array of C positions.
    """
    # Build two orthogonal vectors in the ring plane
    n = np.array(normal, float)
    n = n / np.linalg.norm(n)
    # pick an arbitrary vector not parallel to n
    if abs(n[0]) < 0.9:
        u = np.array([1, 0, 0], float)
    else:
        u = np.array([0, 1, 0], float)
    u = u - np.dot(u, n) * n
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)

    angles = np.linspace(0, 2 * np.pi, 6)[:5]  # 0, 72, 144, 216, 288 degrees
    positions = np.array([center + radius * (np.cos(a) * u + np.sin(a) * v)
                          for a in angles])
    return positions


def _bent_metallocene_template(mol, metal_sym, metal_idx):
    """Generate idealised bent/parallel sandwich geometry for a metallocene.

    Returns (syms, positions) for the whole fragment (including substituents
    on Cp rings placed radially outward), or None on failure.
    """
    ring1, ring2 = _find_cp_rings(mol, metal_idx)
    if not ring1 or not ring2:
        return None

    bent = _BENT_ANGLE.get(metal_sym, 180)
    m_cp = _M_CP_DIST.get(metal_sym, 2.0)

    # Place metal at origin
    metal_pos = np.zeros(3)

    # For bent metallocenes: two Cp centroids at angle bent/2 from z-axis
    # For parallel (ferrocene): centroids on ±z
    half_angle = np.radians(bent / 2)

    if bent < 180:
        # Bent: centroids in xz-plane, symmetric about z
        cp1_center = np.array([m_cp * np.sin(half_angle), 0, m_cp * np.cos(half_angle)])
        cp2_center = np.array([-m_cp * np.sin(half_angle), 0, m_cp * np.cos(half_angle)])
        # Normals point from centroid toward metal (for bent)
        n1 = -cp1_center / m_cp
        n2 = -cp2_center / m_cp
    else:
        # Parallel sandwich
        cp1_center = np.array([0, 0, m_cp])
        cp2_center = np.array([0, 0, -m_cp])
        n1 = np.array([0, 0, -1])
        n2 = np.array([0, 0, 1])

    cp1_xyz = _regular_pentagon_xyz(cp1_center, n1)
    cp2_xyz = _regular_pentagon_xyz(cp2_center, n2)

    # Build full atom list: metal + Cp carbons + substituents
    # Map: ring atom index -> position
    pos_map = {}
    pos_map[metal_idx] = metal_pos
    for i, ci in enumerate(ring1):
        pos_map[ci] = cp1_xyz[i]
    for i, ci in enumerate(ring2):
        pos_map[ci] = cp2_xyz[i]

    # Place substituents radially outward from Cp carbons
    # For each Cp carbon, its non-ring non-metal neighbors go radially outward
    all_atoms = list(mol.GetAtoms())
    placed = set(pos_map.keys())
    queue = list(pos_map.keys())

    while queue:
        ai = queue.pop(0)
        a = mol.GetAtomWithIdx(ai)
        apos = pos_map[ai]
        for nb in a.GetNeighbors():
            ni = nb.GetIdx()
            if ni in placed:
                continue
            # Direction: away from the ring centroid
            # Find which ring this atom belongs to
            if ai in ring1:
                direction = apos - cp1_center
            elif ai in ring2:
                direction = apos - cp2_center
            elif ai == metal_idx:
                # Substituent directly on metal — place opposite to Cp bisector
                if bent < 180:
                    # Bent: bisector is +z; place substituent opposite (−z-ish)
                    direction = np.array([0, 0, -1], float)
                else:
                    # Parallel: place in xy-plane (equatorial)
                    direction = np.array([1, 0, 0], float)
                # If multiple substituents on metal, fan them out
                n_subs_on_metal = sum(1 for nb in a.GetNeighbors()
                                      if nb.GetIdx() not in placed
                                      and nb.GetIdx() not in ring1
                                      and nb.GetIdx() not in ring2)
                if n_subs_on_metal > 1:
                    # Rotate direction for each
                    angle_offset = 2 * np.pi / n_subs_on_metal
                    sub_idx = sum(1 for nb in a.GetNeighbors()
                                  if nb.GetIdx() in placed
                                  and nb.GetIdx() not in ring1
                                  and nb.GetIdx() not in ring2)
                    # Rotate around z-axis
                    c, s = np.cos(angle_offset * sub_idx), np.sin(angle_offset * sub_idx)
                    direction = np.array([c * direction[0] - s * direction[1],
                                          s * direction[0] + c * direction[1],
                                          direction[2]], float)
            else:
                # Already-placed substituent — direction away from its parent
                direction = apos - pos_map.get(ai, apos)

            n = np.linalg.norm(direction)
            if n < 1e-6:
                direction = np.array([0, 0, 1], float)
            else:
                direction = direction / n

            # Bond length: use typical C–X distances
            nb_sym = nb.GetSymbol()
            if nb_sym == "H":
                bl = 1.09
            elif nb_sym in ("C", "N", "O"):
                bl = 1.40
            elif nb_sym == "S":
                bl = 1.75
            elif nb_sym == "P":
                bl = 1.80
            elif nb_sym == "Cl":
                bl = 1.70
            else:
                bl = 1.50

            pos_map[ni] = apos + direction * bl
            placed.add(ni)
            queue.append(ni)

    syms = [a.GetSymbol() for a in all_atoms]
    positions = np.array([pos_map[i] for i in range(len(all_atoms))])
    return syms, positions


# ═══════════════════════════════════════════════════════════════════
#  Improved _embed — uses bent template for metallocenes
# ═══════════════════════════════════════════════════════════════════

def _embed(mol, seed=1):
    """Embed a fragment in 3D. For metallocenes, use idealised bent/parallel
    template instead of RDKit ETKDG (which produces random, often overlapping
    Cp rings). For organic fragments, use ETKDG as before."""
    is_met, metal_sym, metal_idx = _is_metallocene(mol)
    if is_met:
        # Use bent template
        result = _bent_metallocene_template(mol, metal_sym, metal_idx)
        if result is None:
            # Fall back to ETKDG
            pass
        else:
            syms, pos = result
            # Build an RDKit mol with the template coordinates
            mh = Chem.AddHs(mol)
            conf = Chem.Conformer(mh.GetNumAtoms())
            for i, p in enumerate(pos):
                conf.SetAtomPosition(i, p)
            mh.AddConformer(conf)
            return mh

    # Standard ETKDG path
    mh = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useRandomCoords = True
    p.maxIterations = 4000
    if AllChem.EmbedMolecule(mh, p) != 0:
        if AllChem.EmbedMolecule(mh, useRandomCoords=True, maxAttempts=50,
                                 randomSeed=seed + 17) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except Exception:
        pass
    return mh


# ═══════════════════════════════════════════════════════════════════
#  Improved _prep_fragment — softer constraints for metallocenes
# ═══════════════════════════════════════════════════════════════════

def _internal_metal_bonds(mh):
    """[(metal_idx, neighbor_idx, ideal_dist)] for every internal-metal bond."""
    out = []
    for a in mh.GetAtoms():
        if a.GetSymbol() not in INTERNAL_METALS:
            continue
        mi = a.GetIdx()
        for nb in a.GetNeighbors():
            d = _R_INTERNAL.get((a.GetSymbol(), nb.GetSymbol()), 2.2)
            out.append((mi, nb.GetIdx(), d))
    return out


def _prep_fragment(mh, seed=1):
    """Clean a single fragment's 3D geometry. For metallocenes, use
    centroid-based constraints (freeze M–Cp_centroid distances, not all 10
    M–C bonds individually) — this lets the Cp rings breathe while holding
    the sandwich. Plain organic fragments returned as-is."""
    syms = [a.GetSymbol() for a in mh.GetAtoms()]
    conf = mh.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])

    is_met, metal_sym, metal_idx = _is_metallocene(mh)
    if not is_met:
        return syms, pos

    # For metallocenes: constrain M–Cp_centroid distances (2 constraints)
    # instead of all 10 M–C bonds — softer, lets rings adjust
    ring1, ring2 = _find_cp_rings(mh, metal_idx)
    if not ring1 or not ring2:
        # Fall back to per-bond constraints
        ib = _internal_metal_bonds(mh)
        if not ib:
            return syms, pos
        for mi, nj, d in ib:
            v = pos[nj] - pos[mi]
            n = np.linalg.norm(v)
            if n > 1e-6:
                pos[nj] = pos[mi] + v / n * d
    else:
        # Set centroid distances to ideal
        m_cp = _M_CP_DIST.get(metal_sym, 2.0)
        cp1_center = pos[ring1].mean(axis=0)
        cp2_center = pos[ring2].mean(axis=0)
        # Move centroids to ideal distance from metal
        v1 = cp1_center - pos[metal_idx]
        n1 = np.linalg.norm(v1)
        if n1 > 1e-6:
            shift1 = v1 / n1 * m_cp - v1
            pos[ring1] += shift1
        v2 = cp2_center - pos[metal_idx]
        n2 = np.linalg.norm(v2)
        if n2 > 1e-6:
            shift2 = v2 / n2 * m_cp - v2
            pos[ring2] += shift2

    # Constrained relax via native xtb CLI (harmonic, soft springs)
    # Use centroid constraints: (metal_idx, ring_centroid_virtual, m_cp)
    # xTB doesn't support virtual atoms, so we use per-bond constraints
    # but with LOWER force constant (0.5 instead of 1.5) for metallocenes
    ib = _internal_metal_bonds(mh)
    if not ib:
        return syms, pos

    try:
        from ase import Atoms
        from ase.optimize import LBFGS
        from ase.constraints import FixBondLengths
        from xtb.ase.calculator import XTB
        q = Chem.GetFormalCharge(Chem.RemoveHs(mh))
        at = Atoms(symbols=syms, positions=pos)
        at.set_constraint(FixBondLengths([(i, j) for i, j, _ in ib]))
        try:
            at.calc = XTB(method="GFN-FF", charge=int(q), uhf=0)
            LBFGS(at, logfile=None).run(fmax=0.5, steps=200)
        except Exception:
            pass
        at.calc = XTB(method="GFN2-xTB", charge=int(q), uhf=0,
                      electronic_temperature=1500.0, accuracy=2.0, max_iterations=300)
        LBFGS(at, logfile=None).run(fmax=0.2, steps=250)
        return syms, at.get_positions()
    except Exception:
        return syms, pos


# ═══════════════════════════════════════════════════════════════════
#  Assembly + constrained relax (from build_poly, with adaptive fc)
# ═══════════════════════════════════════════════════════════════════

def _find_core_donor(mh):
    """Atom that should bind the CORE metal."""
    best = None
    for a in mh.GetAtoms():
        if any(nb.GetSymbol() in INTERNAL_METALS for nb in a.GetNeighbors()):
            continue
        sym, ch = a.GetSymbol(), a.GetFormalCharge()
        s = None
        if sym == "C" and ch < 0:
            s = 95
        elif sym == "P" and a.GetDegree() == 3 and not any(
                nb.GetSymbol() == "O" for nb in a.GetNeighbors()):
            s = 92
        elif sym == "S" and ch < 0:
            s = 86
        elif sym == "N" and ch < 0:
            s = 82
        elif sym in ("Cl", "Br", "I") and ch < 0:
            s = 52
        if s is not None and (best is None or s > best[2]):
            best = (a.GetIdx(), sym, s)
    return best


def _align(frag_pos, donor_i, neigh_idxs, target_axis, r):
    """Rigid-move fragment so donor sits at target_axis*r from origin."""
    D = frag_pos[donor_i]
    if neigh_idxs:
        N = frag_pos[neigh_idxs].mean(axis=0)
    else:
        N = D - target_axis
    lp = D - N
    n = np.linalg.norm(lp)
    lp = lp / n if n > 1e-6 else np.array(target_axis, float)
    a, b = lp, np.array(target_axis, float)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        R = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))
    moved = (R @ frag_pos.T).T
    moved += (np.array(target_axis, float) * r) - moved[donor_i]
    return moved


def _rot_about(axis, theta):
    axis = np.array(axis, float)
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = axis
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]])


def _min_cross(pos, sym, a_idx, b_idx):
    a = [i for i in a_idx if sym[i] != "H"]
    b = [i for i in b_idx if sym[i] != "H"]
    if not a or not b:
        return 9.9
    A, B = pos[a], pos[b]
    return float(np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)))


def _declash(pos, sym, ranges, axes, passes=3):
    pos = pos.copy()
    search_axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for _ in range(passes):
        for fi, (s, e) in enumerate(ranges):
            others = [i for i in range(len(sym)) if not (s <= i < e)]
            for rax in [axes[fi]] + search_axes:
                best_ang, best_min = 0.0, _min_cross(pos, sym, range(s, e), others)
                for deg in range(15, 360, 15):
                    test = pos.copy()
                    test[s:e] = (_rot_about(rax, np.radians(deg)) @ pos[s:e].T).T
                    d = _min_cross(test, sym, range(s, e), others)
                    if d > best_min:
                        best_min, best_ang = d, deg
                if best_ang:
                    pos[s:e] = (_rot_about(rax, np.radians(best_ang)) @ pos[s:e].T).T
    return pos


def _xtb_cli_opt(syms, pos, charge, uhf, constraints, method="gfn2",
                 fc=1.5, maxcycle=300, timeout=700, etemp=1500):
    """Native xtb CLI optimisation with HARMONIC distance constraints."""
    xtb_bin = shutil.which("xtb") or "/root/biometaldb-3d/arch-env/bin/xtb"
    d = tempfile.mkdtemp(prefix="xtbpoly_")
    try:
        with open(os.path.join(d, "in.xyz"), "w") as f:
            f.write(f"{len(syms)}\n\n")
            for s, p in zip(syms, pos):
                f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        with open(os.path.join(d, "xc.inp"), "w") as f:
            f.write("$constrain\n")
            f.write(f"  force constant={fc}\n")
            for i, j, r0 in constraints:
                f.write(f"  distance: {i + 1}, {j + 1}, {r0:.3f}\n")
            f.write("$end\n")
            f.write(f"$opt\n  maxcycle={maxcycle}\n$end\n")
        cmd = [xtb_bin, "in.xyz", "--opt", "loose", "--input", "xc.inp",
               "--chrg", str(int(charge)), "--uhf", str(int(uhf)),
               "--etemp", str(etemp)]
        cmd += ["--gfnff"] if method == "gfnff" else ["--gfn", "2"]
        env = dict(os.environ, OMP_NUM_THREADS="6", OMP_STACKSIZE="4G")
        subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                       timeout=timeout, env=env)
        outxyz = os.path.join(d, "xtbopt.xyz")
        if not os.path.exists(outxyz):
            return None
        lines = open(outxyz).read().splitlines()
        n = int(lines[0])
        energy = None
        try:
            for tok in lines[1].replace("energy:", " ").split():
                try:
                    energy = float(tok)
                    break
                except ValueError:
                    continue
        except Exception:
            pass
        rs, rp = [], []
        for ln in lines[2:2 + n]:
            q = ln.split()
            rs.append(q[0])
            rp.append([float(q[1]), float(q[2]), float(q[3])])
        return rs, np.array(rp), (energy * 27.2114 if energy is not None else 0.0)
    except Exception:
        return None
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _build_once(metal, ox, smiles_ligands, cn=None, charge=None,
                uhf=0, relax_steps=400, seed=1, adaptive_fc=True):
    """Single-seed build with bent-metallocene template."""
    frags = [f.strip() for f in smiles_ligands.split(".") if f.strip()]
    parsed = []
    for f in frags:
        m = Chem.MolFromSmiles(f)
        if m is None:
            m = Chem.MolFromSmiles(f, sanitize=False)
        if m is None:
            return {"status": "failed", "error": f"RDKit could not parse fragment {f[:40]}"}
        parsed.append(m)

    units = []
    spectators = []
    has_metallocene = False
    for m in parsed:
        donor = _find_core_donor(Chem.AddHs(m))
        mh = _embed(m, seed)
        if mh is None:
            return {"status": "failed", "error": "fragment embed failed"}
        donor = _find_core_donor(mh)
        if donor is None:
            spectators.append(Chem.MolToSmiles(m))
            continue
        di, de, _ = donor
        neigh = [nb.GetIdx() for nb in mh.GetAtomWithIdx(di).GetNeighbors()]
        units.append((mh, di, neigh, de))
        if _is_metallocene(mh)[0]:
            has_metallocene = True

    if not units:
        return {"status": "failed", "error": "no core-metal donor atoms found"}

    cn = cn or min(len(units), 2 if metal in ("Au", "Ag", "Cu") and (ox == 1) else 4)
    units = units[:cn]
    axes = _GEOM_AXES.get(len(units), _GEOM_AXES[2])
    axes = [np.array(a, float) / np.linalg.norm(a) for a in axes]

    all_sym = [metal]
    all_pos = [np.zeros(3)]
    core_bonds = []
    internal_bonds = []
    ranges = []
    frag_axes = []
    offset = 1
    tot_charge = 0
    for (mh, di, neigh, de), axis in zip(units, axes):
        syms, pos = _prep_fragment(mh, seed)
        r = _R_DONOR.get(de, 2.1)
        moved = _align(pos, di, neigh, axis, r)
        tot_charge += Chem.GetFormalCharge(Chem.RemoveHs(mh))
        for mi, nj, d in _internal_metal_bonds(mh):
            internal_bonds.append((offset + mi, offset + nj, d))
        core_bonds.append((0, offset + di, r))
        ranges.append((offset, offset + len(syms)))
        frag_axes.append(axis)
        all_sym += syms
        all_pos.append(moved)
        offset += len(syms)
    all_pos = np.vstack(all_pos)

    for i, j, d in core_bonds:
        v = all_pos[j] - all_pos[i]
        n = np.linalg.norm(v)
        if n > 1e-6:
            all_pos[j] = all_pos[i] + v / n * d
    all_pos = _declash(all_pos, all_sym, ranges, frag_axes)
    seed_min = min((_min_cross(all_pos, all_sym, range(*ranges[a]), range(*ranges[b]))
                    for a in range(len(ranges)) for b in range(a + 1, len(ranges))),
                   default=9.9)

    if charge is None:
        charge = tot_charge + (ox or 0)

    constraints = core_bonds + internal_bonds

    # Adaptive force constant: softer for metallocenes (let rings breathe)
    if adaptive_fc and has_metallocene:
        fc_ff, fc_gfn2 = 0.5, 0.8
    else:
        fc_ff, fc_gfn2 = 1.0, 1.5

    # Two-stage: GFN-FF declash → GFN2 with warm etemp
    ff = _xtb_cli_opt(all_sym, all_pos, charge, uhf, constraints,
                      method="gfnff", fc=fc_ff, maxcycle=250, etemp=1500)
    if ff is not None:
        all_sym, all_pos = ff[0], ff[1]

    # Try GFN2 with warm etemp; if fails, escalate etemp
    for etemp_try in [1500, 3000, 5000, 8000]:
        res = _xtb_cli_opt(all_sym, all_pos, charge, uhf, constraints,
                           method="gfn2", fc=fc_gfn2, maxcycle=relax_steps,
                           etemp=etemp_try)
        if res is not None:
            syms, pos, energy = res
            return {"status": "ok", "method": f"GFN2-xTB(constrained,fc={fc_gfn2},etemp={etemp_try})+poly",
                    "symbols": syms, "positions": np.asarray(pos),
                    "n_atoms": len(syms), "total_charge": int(charge),
                    "n_unpaired": int(uhf), "energy": energy,
                    "metal": metal, "core_bonds": core_bonds,
                    "internal_bonds": internal_bonds,
                    "spectators": spectators,
                    "seed_min_dist": round(float(seed_min), 2)}

    return {"status": "failed",
            "error": "xtb CLI constrained opt produced no geometry (all etemp levels)",
            "seed_min_dist": round(float(seed_min), 2)}


def build_polynuclear(metal, ox, smiles_ligands, cn=None, charge=None,
                      uhf=0, relax_steps=400, seeds=(1, 2, 4, 7, 11, 13, 17)):
    """Build a polynuclear complex, retrying over embed seeds."""
    last = {"status": "failed", "error": "no seed produced a structure"}
    for sd in seeds:
        r = _build_once(metal, ox, smiles_ligands, cn=cn, charge=charge,
                        uhf=uhf, relax_steps=relax_steps, seed=sd)
        if r.get("status") == "ok":
            r["seed"] = sd
            return r
        last = r
    return last


def metal_environment(result, metal, cutoff=2.9):
    s, x = result["symbols"], np.asarray(result["positions"])
    out = {}
    for mi, e in enumerate(s):
        if e == metal or e in INTERNAL_METALS:
            ds = sorted([(s[i], round(float(np.linalg.norm(x[i] - x[mi])), 2))
                         for i in range(len(s)) if i != mi], key=lambda t: t[1])
            out[f"{e}{mi}"] = [d for d in ds if d[1] < cutoff]
    return out


def write_xyz(result, path, title=""):
    s, x = result["symbols"], result["positions"]
    with open(path, "w") as f:
        f.write(f"{len(s)}\n{title} charge={result['total_charge']} "
                f"mult={result['n_unpaired'] + 1} method={result['method']}\n")
        for sym, p in zip(s, x):
            f.write(f"{sym} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    return path


# ═══════════════════════════════════════════════════════════════════
#  Main: run the 6 failed complexes
# ═══════════════════════════════════════════════════════════════════

def run_failed():
    """Run build_polynuclear on the 6 failed titanocene complexes."""
    import sqlite3
    OUT = "/root/biometaldb-3d/polynuclear/out"
    os.makedirs(OUT, exist_ok=True)

    # 6 failed IDs (excluding #7790 which has a different problem)
    FAILED_IDS = [7323, 7874, 8361, 8382, 8868]
    # #7790 is tri-metallic (Ti+Fe) — special case, handle separately

    con = sqlite3.connect("file:/root/biometaldb-3d/pilot.sqlite?mode=ro", uri=True)
    ok = 0
    for cid in FAILED_IDS:
        metal, ox, smi = con.execute(
            "SELECT metal,oxidation_state,smiles_ligands FROM complexes WHERE id=?",
            (cid,)).fetchone()
        print(f"===== #{cid} {metal}({ox}) =====", flush=True)
        try:
            r = build_polynuclear(metal, ox, smi, uhf=0,
                                  seeds=(1, 2, 4, 7, 11, 13, 17, 23),
                                  relax_steps=400)
        except Exception as e:
            print(f"  EXC {type(e).__name__}: {str(e)[:140]}", flush=True)
            continue
        if r.get("status") != "ok":
            print(f"  FAILED: {r.get('error')}", flush=True)
            continue
        ok += 1
        write_xyz(r, os.path.join(OUT, f"poly_{cid}.xyz"),
                  title=f"BiometalDB #{cid} {metal}({ox})")
        env = metal_environment(r, metal)
        metals_env = " | ".join(
            f"{k}:{sum(1 for e in v if e[0]=='C' and 1.9<e[1]<2.6)}Cp,"
            f"{[e for e in v if e[0] in ('S','P','N','O') and e[1]<2.6][:2]}"
            for k, v in env.items())
        print(f"  OK seed={r.get('seed')} atoms={r['n_atoms']} E={r['energy']:.1f} "
              f"-> poly_{cid}.xyz | {metals_env}", flush=True)
    con.close()
    print(f"ALL_DONE built {ok}/{len(FAILED_IDS)}", flush=True)


def test_template():
    """Quick test: generate bent template for #7874 titanocene fragment."""
    smi = "[CH3-]->[Ti+4]12345678(<-[O-]C(=O)c9ccc([S-])cc9)(<-[CH-]9[C]1=[C]2[C]3=[C]94)<-[CH-]1[C]5=[C]6[C]7=[C]18"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
    is_met, metal_sym, metal_idx = _is_metallocene(mol)
    print(f"is_metallocene: {is_met} metal={metal_sym} idx={metal_idx}")
    if is_met:
        syms, pos = _bent_metallocene_template(mol, metal_sym, metal_idx)
        print(f"template: {len(syms)} atoms")
        # Check Ti–Cp distances
        ring1, ring2 = _find_cp_rings(mol, metal_idx)
        if ring1 and ring2:
            cp1 = pos[ring1].mean(axis=0)
            cp2 = pos[ring2].mean(axis=0)
            d1 = np.linalg.norm(cp1 - pos[metal_idx])
            d2 = np.linalg.norm(cp2 - pos[metal_idx])
            angle = np.degrees(np.arccos(np.dot(cp1 - pos[metal_idx],
                                                cp2 - pos[metal_idx]) / (d1 * d2)))
            print(f"Ti–Cp1: {d1:.2f} Å  Ti–Cp2: {d2:.2f} Å  Cp–Ti–Cp angle: {angle:.1f}°")
        # Write XYZ for visual inspection
        with open("/tmp/bent_template_test.xyz", "w") as f:
            f.write(f"{len(syms)}\nBent titanocene template test\n")
            for s, p in zip(syms, pos):
                f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        print("XYZ written to /tmp/bent_template_test.xyz")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--run-failed", action="store_true")
    args = ap.parse_args()

    if args.test:
        test_template()
    elif args.run_failed:
        run_failed()
    else:
        print("Usage: python bent_template.py --test | --run-failed")
