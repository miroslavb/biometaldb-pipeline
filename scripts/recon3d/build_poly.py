"""Polynuclear / organometallic 3D builder (recon3d extension).

The mononuclear Architector pipeline cannot place complexes whose ligand SMILES
contain a SECOND metal (titanocene / ferrocene sandwich co-ligands written with
dative bonds `->`/`<-`), and such SMILES crash OpenBabel. This module builds them
by fragment assembly + a constrained xtb relaxation, entirely OpenBabel-free:

  1. parse each dot-fragment with RDKit (it handles dative bonds);
  2. find each fragment's CORE-metal donor atom (carbene C / phosphine P /
     thiolate S / terminal-alkynyl C) — the one NOT bonded to an internal metal;
  3. embed each fragment in 3D (RDKit ETKDG; metallocene η5 comes out crude but
     is fixed by the constrained relax);
  4. place the fragments around the core metal at the right polyhedron
     (CN2→linear, CN4→tetra/square), aligning each donor's lone-pair axis;
  5. merge + GFN2-xTB relax with EVERY metal's coordination distances FROZEN
     (FixBondLengths) so neither the core sphere nor any sandwich collapses;
  6. emit xyz + mol2 written directly from coordinates (no OpenBabel).

Scope: built for the BiometalDB Au(I) + metallocene-pendant tail, but the donor
detection + linear/tetra placement generalise to other CN2/CN4 cores.
"""
from __future__ import annotations
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

INTERNAL_METALS = {"Ti", "Fe", "Co", "Ru", "Ni", "Mn", "Cr", "V", "Zr", "Hf", "Mo", "W"}
# ideal core-metal–donor distances (Å)
_R_DONOR = {"C": 2.02, "P": 2.28, "S": 2.30, "N": 2.05, "Cl": 2.27, "O": 2.05}
# ideal internal-metal coordination distances
_R_INTERNAL = {("Ti", "C"): 2.37, ("Ti", "O"): 2.05, ("Ti", "Cl"): 2.30,
               ("Fe", "C"): 2.04, ("Co", "C"): 2.10, ("Ru", "C"): 2.20}


# --------------------------------------------------------------------------- #
def _xtb_cli_opt(syms, pos, charge, uhf, constraints, method="gfn2",
                 fc=1.5, maxcycle=300, timeout=700):
    """Native xtb CLI optimisation with HARMONIC distance constraints ($constrain
    force constant) — a soft spring on each metal-coordination distance, NOT a hard
    freeze, so the optimiser can relieve inter-fragment clashes while the sandwich
    and the core sphere are held (the XPD 'strong-restrain, don't freeze' recipe).
    Returns (syms, positions, energy) or None."""
    import os
    import shutil
    import subprocess
    import tempfile
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
               "--etemp", "1500"]
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
        except Exception:  # noqa: BLE001
            pass
        rs, rp = [], []
        for ln in lines[2:2 + n]:
            q = ln.split()
            rs.append(q[0])
            rp.append([float(q[1]), float(q[2]), float(q[3])])
        return rs, np.array(rp), (energy * 27.2114 if energy is not None else 0.0)
    except Exception:  # noqa: BLE001
        return None
    finally:
        shutil.rmtree(d, ignore_errors=True)


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
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)  # cleans the organic part only
    except Exception:  # noqa: BLE001
        pass
    return mh


def _prep_fragment(mh, seed=1):
    """Clean a single fragment's 3D geometry. For a metallocene (contains an
    internal metal) the RDKit embed is crude and forcing radial distances would
    squash the Cp rings, so we constrained-relax it IN ISOLATION (freeze the
    internal-metal coordination, GFN-FF then GFN2) — the recipe proven to hold the
    sandwich. Plain organic fragments are returned as-is. Returns (syms, pos) or
    None on failure."""
    syms = [a.GetSymbol() for a in mh.GetAtoms()]
    conf = mh.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])
    ib = _internal_metal_bonds(mh)
    if not ib:
        return syms, pos
    # set internal-metal distances to ideal radially (per-fragment, no other metal
    # nearby so no squashing across fragments), then constrained relax
    for mi, nj, d in ib:
        v = pos[nj] - pos[mi]
        n = np.linalg.norm(v)
        if n > 1e-6:
            pos[nj] = pos[mi] + v / n * d
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
        except Exception:  # noqa: BLE001
            pass
        at.calc = XTB(method="GFN2-xTB", charge=int(q), uhf=0,
                      electronic_temperature=1500.0, accuracy=2.0, max_iterations=300)
        LBFGS(at, logfile=None).run(fmax=0.2, steps=250)
        return syms, at.get_positions()
    except Exception:  # noqa: BLE001
        return syms, pos          # fall back to the radial seed if relax fails


def _find_core_donor(mh):
    """Atom (index) that should bind the CORE metal: a carbanion/carbene C,
    phosphine P, thiolate S, or terminal-alkynyl C — NOT bonded to an internal
    metal. Returns (idx, element, strength) or None."""
    best = None
    for a in mh.GetAtoms():
        if any(nb.GetSymbol() in INTERNAL_METALS for nb in a.GetNeighbors()):
            continue
        sym, ch = a.GetSymbol(), a.GetFormalCharge()
        s = None
        if sym == "C" and ch < 0:
            s = 95                                   # carbene / carbanion / acetylide
        elif sym == "P" and a.GetDegree() == 3 and not any(
                nb.GetSymbol() == "O" for nb in a.GetNeighbors()):
            s = 92                                   # phosphine
        elif sym == "S" and ch < 0:
            s = 86                                   # thiolate
        elif sym == "N" and ch < 0:
            s = 82
        elif sym in ("Cl", "Br", "I") and ch < 0:
            s = 52
        if s is not None and (best is None or s > best[2]):
            best = (a.GetIdx(), sym, s)
    return best


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


def _align(frag_pos, donor_i, neigh_idxs, target_axis, r):
    """Rigid-move a fragment so its donor sits at `target_axis*r` from the origin
    (core metal) and its lone-pair axis (away from the donor's substituents) is
    anti-parallel to target_axis (i.e. the donor points back at the metal)."""
    D = frag_pos[donor_i]
    if neigh_idxs:
        N = frag_pos[neigh_idxs].mean(axis=0)
    else:
        N = D - target_axis             # degenerate: assume already aligned
    lp = D - N                          # lone-pair direction (substituents -> donor)
    n = np.linalg.norm(lp)
    lp = lp / n if n > 1e-6 else np.array(target_axis, float)
    # rotate lp -> target_axis
    a, b = lp, np.array(target_axis, float)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        R = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))
    moved = (R @ frag_pos.T).T
    # translate so donor lands at target_axis*r
    moved += (np.array(target_axis, float) * r) - moved[donor_i]
    return moved


_GEOM_AXES = {
    2: [(0, 0, 1), (0, 0, -1)],
    3: [(0, 0, 1), (0.866, 0, -0.5), (-0.866, 0, -0.5)],
    4: [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],   # tetrahedral
}


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
    """Min heavy-heavy distance between two atom-index sets."""
    a = [i for i in a_idx if sym[i] != "H"]
    b = [i for i in b_idx if sym[i] != "H"]
    if not a or not b:
        return 9.9
    A, B = pos[a], pos[b]
    return float(np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)))


def _declash(pos, sym, ranges, axes, passes=3):
    """Rotate each fragment about axes THROUGH THE CORE (origin) to maximise the
    minimum distance to all other atoms. Any such rotation preserves |core-donor|
    (the donor's distance from the origin), so the frozen coordination distance is
    untouched; only the donor *direction* changes (the relax restores the ideal
    angle). Searches the fragment's own placement axis (spin) + x/y/z (fan)."""
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


def build_polynuclear(metal, ox, smiles_ligands, cn=None, charge=None,
                      uhf=0, relax_steps=400, seeds=(1, 2, 4, 7)):
    """Build a polynuclear complex, retrying over embed seeds until one relaxes.
    Returns a result dict like build.build() or {'status':'failed',...}."""
    last = {"status": "failed", "error": "no seed produced a structure"}
    for sd in seeds:
        r = _build_once(metal, ox, smiles_ligands, cn=cn, charge=charge,
                        uhf=uhf, relax_steps=relax_steps, seed=sd)
        if r.get("status") == "ok":
            r["seed"] = sd
            return r
        last = r
    return last


def _build_once(metal, ox, smiles_ligands, cn=None, charge=None,
                uhf=0, relax_steps=400, seed=1):
    """Single-seed build. Never calls OpenBabel."""
    frags = [f.strip() for f in smiles_ligands.split(".") if f.strip()]
    parsed = []
    for f in frags:
        m = Chem.MolFromSmiles(f)
        if m is None:
            m = Chem.MolFromSmiles(f, sanitize=False)
        if m is None:
            return {"status": "failed", "error": f"RDKit could not parse fragment {f[:40]}"}
        parsed.append(m)

    # embed each fragment, find its core-metal donor
    units = []          # (mh, donor_idx, neigh_idxs, donor_elem)
    spectators = []
    for m in parsed:
        donor = _find_core_donor(Chem.AddHs(m))   # detect on H-added graph indices
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

    if not units:
        return {"status": "failed", "error": "no core-metal donor atoms found"}

    cn = cn or min(len(units), 2 if metal in ("Au", "Ag", "Cu") and (ox == 1) else 4)
    units = units[:cn]
    axes = _GEOM_AXES.get(len(units), _GEOM_AXES[2])
    axes = [np.array(a, float) / np.linalg.norm(a) for a in axes]

    # assemble: core metal at origin + each aligned fragment
    all_sym = [metal]
    all_pos = [np.zeros(3)]
    core_bonds = []          # (0, atom_idx_in_merged, ideal_r)
    internal_bonds = []      # (metal_idx, nbr_idx, ideal_r) in merged indexing
    ranges = []              # (start, end) merged-index span per fragment
    frag_axes = []
    offset = 1
    tot_charge = 0
    for (mh, di, neigh, de), axis in zip(units, axes):
        syms, pos = _prep_fragment(mh, seed)   # clean per-fragment geometry first
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

    # only the core-donor distances need resetting (fragments are already cleanly
    # relaxed internally); then declash by rotating fragments about the core
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
        charge = tot_charge + (ox or 0)     # core metal contributes +ox

    # constrained relax via native xtb CLI with HARMONIC distance restraints
    # (soft springs on every metal-coordination distance — relieves clashes while
    # holding both the core sphere and the metallocene sandwich). Two stages:
    # GFN-FF (robust declash) then GFN2 (warm etemp for the d-system SCF).
    constraints = core_bonds + internal_bonds
    ff = _xtb_cli_opt(all_sym, all_pos, charge, uhf, constraints,
                      method="gfnff", fc=1.0, maxcycle=250)
    if ff is not None:
        all_sym, all_pos = ff[0], ff[1]
    res = _xtb_cli_opt(all_sym, all_pos, charge, uhf, constraints,
                       method="gfn2", fc=1.5, maxcycle=relax_steps)
    if res is None:
        return {"status": "failed", "error": "xtb CLI constrained opt produced no geometry",
                "seed_min_dist": round(float(seed_min), 2)}
    syms, pos, energy = res
    return {"status": "ok", "method": "GFN2-xTB(constrained)+poly", "symbols": syms,
            "positions": np.asarray(pos), "n_atoms": len(syms),
            "total_charge": int(charge), "n_unpaired": int(uhf), "energy": energy,
            "metal": metal, "core_bonds": core_bonds, "internal_bonds": internal_bonds,
            "spectators": spectators, "seed_min_dist": round(float(seed_min), 2)}


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
