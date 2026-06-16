"""Fixed bent-metallocene template: use RDKit ETKDG for the whole fragment,
then replace only the Cp-ring carbon positions with idealised bent/parallel
geometry.  Substituents stay where RDKit placed them (which is chemically
reasonable for organic groups), only the sandwich core is corrected.

This avoids the zero-distance collapse of the previous pure-template approach
while still fixing the root cause (ETKDG's random Cp-ring placement)."""

import warnings; warnings.filterwarnings("ignore")
import sys, os, subprocess, tempfile, time, argparse, json, shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

_BENT_ANGLE = {"Ti": 130, "Zr": 135, "Hf": 135}
_SANDWICH_METALS = {"Ti", "Zr", "Hf", "Fe", "Co", "Ru", "Cr", "Mn", "Ni", "V"}
_M_CP_DIST = {"Ti": 2.05, "Zr": 2.20, "Hf": 2.18,
              "Fe": 1.65, "Co": 1.70, "Ru": 1.80,
              "Cr": 1.70, "Mn": 1.75, "Ni": 1.70, "V": 1.90}
_CP_RADIUS = 1.21

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


def _is_metallocene(mol):
    for a in mol.GetAtoms():
        if a.GetSymbol() not in _SANDWICH_METALS:
            continue
        c_neighbors = [nb for nb in a.GetNeighbors() if nb.GetSymbol() == "C"]
        if len(c_neighbors) >= 8:
            return True, a.GetSymbol(), a.GetIdx()
    return False, None, None


def _find_cp_rings(mol, metal_idx):
    metal = mol.GetAtomWithIdx(metal_idx)
    c_neighbors = [nb.GetIdx() for nb in metal.GetNeighbors() if nb.GetSymbol() == "C"]
    if len(c_neighbors) < 8:
        return [], []
    adj = {i: set() for i in c_neighbors}
    for i in c_neighbors:
        ai = mol.GetAtomWithIdx(i)
        for nb in ai.GetNeighbors():
            if nb.GetIdx() in c_neighbors:
                adj[i].add(nb.GetIdx())
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


def _regular_pentagon_xyz(center, normal, radius=_CP_RADIUS):
    n = np.array(normal, float)
    n = n / np.linalg.norm(n)
    if abs(n[0]) < 0.9:
        u = np.array([1, 0, 0], float)
    else:
        u = np.array([0, 1, 0], float)
    u = u - np.dot(u, n) * n
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    angles = np.linspace(0, 2 * np.pi, 6)[:5]
    return np.array([center + radius * (np.cos(a) * u + np.sin(a) * v)
                     for a in angles])


def _correct_sandwich_geometry(mh, metal_sym, metal_idx):
    """Replace Cp-ring carbon positions with idealised bent/parallel geometry.
    All other atoms (substituents, H, etc.) stay where RDKit placed them,
    then are rigidly moved along with their parent Cp ring.

    Returns (syms, positions) or None."""
    ring1, ring2 = _find_cp_rings(mh, metal_idx)
    if not ring1 or not ring2:
        return None

    bent = _BENT_ANGLE.get(metal_sym, 180)
    m_cp = _M_CP_DIST.get(metal_sym, 2.0)
    half_angle = np.radians(bent / 2)

    # Current positions
    conf = mh.GetConformer()
    old_pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])
    old_metal = old_pos[metal_idx].copy()

    # Current Cp centroids
    old_cp1 = old_pos[ring1].mean(axis=0)
    old_cp2 = old_pos[ring2].mean(axis=0)

    # Ideal Cp centroids
    if bent < 180:
        new_cp1 = np.array([m_cp * np.sin(half_angle), 0, m_cp * np.cos(half_angle)])
        new_cp2 = np.array([-m_cp * np.sin(half_angle), 0, m_cp * np.cos(half_angle)])
        n1 = -new_cp1 / m_cp
        n2 = -new_cp2 / m_cp
    else:
        new_cp1 = np.array([0, 0, m_cp])
        new_cp2 = np.array([0, 0, -m_cp])
        n1 = np.array([0, 0, -1])
        n2 = np.array([0, 0, 1])

    # Ideal Cp carbon positions
    new_cp1_xyz = _regular_pentagon_xyz(new_cp1, n1)
    new_cp2_xyz = _regular_pentagon_xyz(new_cp2, n2)

    # Build new positions: start with old, then move Cp rings + their
    # substituents rigidly to the ideal geometry
    new_pos = old_pos.copy()

    # Move metal to origin
    new_pos -= old_metal

    # For each ring: find the rigid transform that maps old ring centroid
    # → new ring centroid, and old ring plane → new ring plane.
    # Then apply to all atoms in that ring's "subtree" (ring carbons +
    # everything attached to them that isn't the metal or the other ring).

    def find_subtree(ring_indices, metal_idx, other_ring):
        """All atoms attached to this ring (transitively), excluding metal and other ring."""
        subtree = set(ring_indices)
        queue = list(ring_indices)
        while queue:
            ai = queue.pop(0)
            a = mh.GetAtomWithIdx(ai)
            for nb in a.GetNeighbors():
                ni = nb.GetIdx()
                if ni == metal_idx or ni in other_ring or ni in subtree:
                    continue
                subtree.add(ni)
                queue.append(ni)
        return subtree

    # Ring 1
    sub1 = find_subtree(ring1, metal_idx, set(ring2))
    # Compute rigid transform: old_cp1_center -> new_cp1, with ring plane alignment
    old_cp1_c = old_pos[list(ring1)].mean(axis=0) - old_metal
    # Translation
    t1 = new_cp1 - old_cp1_c
    # Rotation: align old ring plane normal to new
    # Old normal: cross product of two in-plane vectors
    old_v1 = old_pos[ring1[1]] - old_pos[ring1[0]]
    old_v2 = old_pos[ring1[2]] - old_pos[ring1[0]]
    old_n1 = np.cross(old_v1, old_v2)
    old_n1 = old_n1 / np.linalg.norm(old_n1)
    # Rotation from old_n1 to n1
    a, b = old_n1, n1
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        R1 = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R1 = np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))
    # Apply to subtree: translate to origin, rotate, translate to new position
    for i in sub1:
        new_pos[i] = R1 @ (old_pos[i] - old_metal - old_cp1_c) + new_cp1

    # Now place Cp carbons at exact ideal positions
    for i, ci in enumerate(ring1):
        new_pos[ci] = new_cp1_xyz[i]

    # Ring 2
    sub2 = find_subtree(ring2, metal_idx, set(ring1))
    old_cp2_c = old_pos[list(ring2)].mean(axis=0) - old_metal
    t2 = new_cp2 - old_cp2_c
    old_v1 = old_pos[ring2[1]] - old_pos[ring2[0]]
    old_v2 = old_pos[ring2[2]] - old_pos[ring2[0]]
    old_n2 = np.cross(old_v1, old_v2)
    old_n2 = old_n2 / np.linalg.norm(old_n2)
    a, b = old_n2, n2
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        R2 = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R2 = np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))
    for i in sub2:
        new_pos[i] = R2 @ (old_pos[i] - old_metal - old_cp2_c) + new_cp2
    for i, ci in enumerate(ring2):
        new_pos[ci] = new_cp2_xyz[i]

    # Metal stays at origin
    new_pos[metal_idx] = np.zeros(3)

    syms = [a.GetSymbol() for a in mh.GetAtoms()]
    return syms, new_pos


# ═══════════════════════════════════════════════════════════════════
#  Improved _embed — ETKDG + sandwich correction
# ═══════════════════════════════════════════════════════════════════

def _embed(mol, seed=1):
    """Embed via RDKit ETKDG, then correct metallocene sandwich geometry.
    Dative bonds (->, <-) are already converted to single bonds by _build_once,
    so RDKit can embed all fragments normally."""
    mh = Chem.AddHs(mol)

    # Standard ETKDG path for all fragments
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useRandomCoords = True
    p.maxIterations = 4000
    if AllChem.EmbedMolecule(mh, p) != 0:
        if AllChem.EmbedMolecule(mh, useRandomCoords=True, maxAttempts=200,
                                 randomSeed=seed + 17) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except Exception:
        pass

    # Correct sandwich if metallocene
    is_met, metal_sym, metal_idx = _is_metallocene(mh)
    if is_met:
        result = _correct_sandwich_geometry(mh, metal_sym, metal_idx)
        if result is not None:
            syms, pos = result
            conf = Chem.Conformer(mh.GetNumAtoms())
            for i, p in enumerate(pos):
                conf.SetAtomPosition(i, p)
            mh.RemoveAllConformers()
            mh.AddConformer(conf)

    return mh


# ═══════════════════════════════════════════════════════════════════
#  Fragment prep + assembly (same as before, with adaptive fc)
# ═══════════════════════════════════════════════════════════════════

def _internal_metal_bonds(mh):
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
    syms = [a.GetSymbol() for a in mh.GetAtoms()]
    conf = mh.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])
    is_met, metal_sym, metal_idx = _is_metallocene(mh)
    if not is_met:
        return syms, pos

    # For metallocenes: constrained relax with centroid constraints
    ring1, ring2 = _find_cp_rings(mh, metal_idx)
    ib = _internal_metal_bonds(mh)
    if not ib:
        return syms, pos

    # Set centroid distances to ideal
    if ring1 and ring2:
        m_cp = _M_CP_DIST.get(metal_sym, 2.0)
        cp1_center = pos[ring1].mean(axis=0)
        cp2_center = pos[ring2].mean(axis=0)
        v1 = cp1_center - pos[metal_idx]
        n1 = np.linalg.norm(v1)
        if n1 > 1e-6:
            pos[ring1] += v1 / n1 * m_cp - v1
        v2 = cp2_center - pos[metal_idx]
        n2 = np.linalg.norm(v2)
        if n2 > 1e-6:
            pos[ring2] += v2 / n2 * m_cp - v2

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


def _find_core_donor(mh):
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
    frags = [f.strip() for f in smiles_ligands.split(".") if f.strip()]
    parsed = []
    for f in frags:
        # Convert dative bonds (->, <-) to single bonds so RDKit can
        # embed the fragment.  The sandwich geometry is corrected later
        # by _correct_sandwich_geometry; the bond type doesn't matter
        # for the xTB constrained relax.
        f_fixed = f.replace("->", "-").replace("<-", "-")
        m = Chem.MolFromSmiles(f_fixed)
        if m is None:
            m = Chem.MolFromSmiles(f_fixed, sanitize=False)
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

    if adaptive_fc and has_metallocene:
        fc_ff, fc_gfn2 = 0.5, 0.8
    else:
        fc_ff, fc_gfn2 = 1.0, 1.5

    ff = _xtb_cli_opt(all_sym, all_pos, charge, uhf, constraints,
                      method="gfnff", fc=fc_ff, maxcycle=250, etemp=1500)
    if ff is not None:
        all_sym, all_pos = ff[0], ff[1]

    for etemp_try in [1500, 3000, 5000, 8000]:
        res = _xtb_cli_opt(all_sym, all_pos, charge, uhf, constraints,
                           method="gfn2", fc=fc_gfn2, maxcycle=relax_steps,
                           etemp=etemp_try)
        if res is not None:
            syms, pos, energy = res
            return {"status": "ok",
                    "method": f"GFN2-xTB(constrained,fc={fc_gfn2},etemp={etemp_try})+poly",
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
#  Main
# ═══════════════════════════════════════════════════════════════════

def run_failed():
    import sqlite3
    OUT = "/root/biometaldb-3d/polynuclear/out"
    os.makedirs(OUT, exist_ok=True)
    FAILED_IDS = [7323, 7874, 8361, 8382, 8868]
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
    """Quick test: generate corrected sandwich for #7874 titanocene fragment."""
    smi = "[CH3-]->[Ti+4]12345678(<-[O-]C(=O)c9ccc([S-])cc9)(<-[CH-]9[C]1=[C]2[C]3=[C]94)<-[CH-]1[C]5=[C]6[C]7=[C]18"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
    mh = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = 1
    p.useRandomCoords = True
    p.maxIterations = 4000
    AllChem.EmbedMolecule(mh, p)
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except:
        pass

    is_met, metal_sym, metal_idx = _is_metallocene(mh)
    print(f"is_metallocene: {is_met} metal={metal_sym} idx={metal_idx}")
    if is_met:
        syms, pos = _correct_sandwich_geometry(mh, metal_sym, metal_idx)
        print(f"corrected: {len(syms)} atoms")
        ring1, ring2 = _find_cp_rings(mh, metal_idx)
        if ring1 and ring2:
            cp1 = pos[ring1].mean(axis=0)
            cp2 = pos[ring2].mean(axis=0)
            d1 = np.linalg.norm(cp1 - pos[metal_idx])
            d2 = np.linalg.norm(cp2 - pos[metal_idx])
            angle = np.degrees(np.arccos(
                np.dot(cp1 - pos[metal_idx], cp2 - pos[metal_idx]) / (d1 * d2)))
            print(f"Ti–Cp1: {d1:.2f} Å  Ti–Cp2: {d2:.2f} Å  Cp–Ti–Cp angle: {angle:.1f}°")
        # Check for zero-distance overlaps
        from itertools import combinations
        min_d = 999
        for i, j in combinations(range(len(syms)), 2):
            d = np.linalg.norm(pos[i] - pos[j])
            if d < 0.01:
                print(f"  ZERO-DISTANCE: {syms[i]}{i}–{syms[j]}{j} = {d:.4f} Å")
            if d < min_d:
                min_d = d
        print(f"Min interatomic distance: {min_d:.4f} Å")
        with open("/tmp/bent_template_v2_test.xyz", "w") as f:
            f.write(f"{len(syms)}\nCorrected sandwich test\n")
            for s, p in zip(syms, pos):
                f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        print("XYZ → /tmp/bent_template_v2_test.xyz")


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
        print("Usage: python bent_template_v2.py --test | --run-failed")
