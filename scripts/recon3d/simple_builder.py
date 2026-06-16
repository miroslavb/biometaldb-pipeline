#!/usr/bin/env python3
"""Simple 3D builder for complexes that fail in Architector/xTB pipeline.

Strategy: NO xTB constrained relax — just RDKit 3D embedding + custom
metal-coordination placement + MMFF cleanup. Much faster and more robust
for "hard" complexes (CN=2/4/6 with monodentate or bidentate ligands).

Covers:
- haptic_overcrowded: Au/Ag/Cu get haptic-filtered (only Fe/Ru/Os allowed)
- insufficient_donors: auto-add [Cl-]/[OH-]/[NH3] to reach CN for Au(I/III)
- corrupt_source_donors: parse in sanitize=False
- geometry_unassemblable: simple metal-centered assembly without xTB
"""
import warnings; warnings.filterwarnings('ignore')
import sys, os, argparse, json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# Donor radii (Å) by element
_R_DONOR = {"C": 2.02, "P": 2.28, "S": 2.30, "N": 2.05, "Cl": 2.27, "O": 2.05}

# Haptic metals (Fe, Ru, Os do η5/η6; Au/Ag/Cu do not)
HAPTIC_METALS = {"Fe", "Ru", "Os", "Co", "Rh", "Ir", "Mn", "Cr", "V"}
HAPTIC_ATOMS = {"C", "N"}   # haptic hapticity determined by M-C/N count

# Axial layouts
_LINEAR = [(0, 0, 1), (0, 0, -1)]
_SQ_PLANAR = [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]
_TETRA = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
_OCTA = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

# CN-by-metal-ox expected coordination number
EXPECTED_CN = {
    ("Au", 1): 2, ("Au", 3): 4,
    ("Ag", 1): 2, ("Ag", 3): 4,
    ("Cu", 1): 2, ("Cu", 2): 4, ("Cu", 3): 4,
    ("Pt", 2): 4, ("Pt", 4): 6,
    ("Pd", 2): 4, ("Pd", 4): 6,
    ("Rh", 1): 4, ("Rh", 3): 6,
    ("Ir", 1): 4, ("Ir", 3): 6,
    ("Ru", 2): 6, ("Ru", 3): 6,
    ("Os", 2): 6, ("Os", 3): 6,
    ("Co", 2): 6, ("Co", 3): 6,
    ("Fe", 2): 6, ("Fe", 3): 6,
    ("Ni", 2): 4, ("Ni", 0): 4,
    ("Zn", 2): 4,
    ("Mn", 2): 6,
}


def _is_haptic_fragment(mol, metal_sym):
    """True if fragment contains a haptic (η5/η6) metal-carbon sandwich.
    Returns (is_haptic, metal_idx) or (False, None)."""
    if metal_sym not in HAPTIC_METALS:
        return False, None
    for a in mol.GetAtoms():
        if a.GetSymbol() != metal_sym: continue
        c_neighbors = [nb for nb in a.GetNeighbors() if nb.GetSymbol() == "C"]
        # η5-Cp = 5 C from one ring; η6-arene = 6 C from one ring
        if len(c_neighbors) >= 5:
            return True, a.GetIdx()
    return False, None


def _find_donor_simple(mh, metal_sym):
    """Find best donor atom for the metal.
    Priority: carbanion C > phosphine P > thiolate S > neutral N (bidentate) > halide.
    """
    best = None
    for a in mh.GetAtoms():
        sym, ch = a.GetSymbol(), a.GetFormalCharge()
        s = None
        if sym == "C" and ch < 0: s = 95
        elif sym == "P" and a.GetDegree() == 3 and not any(
                nb.GetSymbol() == "O" for nb in a.GetNeighbors()): s = 92
        elif sym == "S" and ch < 0: s = 86
        elif sym == "N" and ch < 0: s = 82
        # bidentate N (e.g. phenanthroline, bipy) — neutral but chelating
        elif sym == "N" and a.GetDegree() == 2 and a.GetIsAromatic():
            s = 70
        elif sym in ("Cl", "Br", "I") and ch < 0: s = 52
        elif sym == "O" and a.GetDegree() >= 1 and ch < 0: s = 60
        if s is not None and (best is None or s > best[2]):
            best = (a.GetIdx(), sym, s)
    return best


def _align(frag_pos, donor_i, neigh_idxs, target_axis, r):
    D = frag_pos[donor_i]
    N = frag_pos[neigh_idxs].mean(axis=0) if neigh_idxs else D - target_axis
    lp = D - N
    lp = lp / np.linalg.norm(lp) if np.linalg.norm(lp) > 1e-6 else np.array(target_axis, float)
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


def _min_cross(pos, sym, a_idx, b_idx):
    a = [i for i in a_idx if sym[i] != "H"]
    b = [i for i in b_idx if sym[i] != "H"]
    if not a or not b: return 9.9
    return float(np.min(np.linalg.norm(pos[a][:, None, :] - pos[b][None, :, :], axis=-1)))


def _declash(pos, sym, ranges, axes, passes=3):
    """Rotate fragments about the origin to relieve inter-fragment clashes."""
    pos = pos.copy()
    search_axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    def _rot(a, t):
        a = np.array(a, float); a = a / np.linalg.norm(a)
        c, s = np.cos(t), np.sin(t)
        x, y, z = a
        return np.array([
            [c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
            [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
            [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])
    for _ in range(passes):
        for fi, (s, e) in enumerate(ranges):
            others = [i for i in range(len(sym)) if not (s <= i < e)]
            for rax in [axes[fi]] + search_axes:
                best_ang, best_min = 0.0, _min_cross(pos, sym, range(s, e), others)
                for deg in range(15, 360, 15):
                    test = pos.copy()
                    test[s:e] = (_rot(rax, np.radians(deg)) @ pos[s:e].T).T
                    d = _min_cross(test, sym, range(s, e), others)
                    if d > best_min: best_min, best_ang = d, deg
                if best_ang: pos[s:e] = (_rot(rax, np.radians(best_ang)) @ pos[s:e].T).T
    return pos


def _embed(mol, seed=1):
    """ETKDG embed. Convert dative bonds (->, <-) to single first.
    Skip MMFF for large mols (slow) — ETKDG geometry is good enough."""
    smi = Chem.MolToSmiles(mol)
    smi = smi.replace("->", "-").replace("<-", "-")
    m = Chem.MolFromSmiles(smi) or Chem.MolFromSmiles(smi, sanitize=False)
    if m is None: return None
    mh = Chem.AddHs(m)
    n_heavy = mh.GetNumHeavyAtoms()
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useRandomCoords = True
    p.maxIterations = 4000
    if AllChem.EmbedMolecule(mh, p) != 0:
        if AllChem.EmbedMolecule(mh, useRandomCoords=True, maxAttempts=200,
                                 randomSeed=seed + 17) != 0:
            return None
    # Skip MMFF for big mols — saves ~5-15s each
    if n_heavy < 100:
        try: AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
        except: pass
    return mh


def build(metal, ox, smiles_ligands, max_atoms=400):
    """Build a 3D structure. Returns dict with 'status' and 'symbols'/'positions' on success."""
    frags = [f.strip() for f in smiles_ligands.split(".") if f.strip()]
    parsed = []
    for f in frags:
        m = Chem.MolFromSmiles(f) or Chem.MolFromSmiles(f, sanitize=False)
        if m is None: continue
        parsed.append(m)

    if not parsed:
        return {"status": "failed", "error": "no fragments parsed"}

    units = []
    spectators = []
    for m in parsed:
        # Skip haptic fragments for non-haptic metals (Au, Ag, Cu)
        is_hap, _ = _is_haptic_fragment(m, metal)
        if is_hap and metal not in HAPTIC_METALS:
            spectators.append(Chem.MolToSmiles(m))
            continue
        mh = _embed(m, 1)
        if mh is None or mh.GetNumConformers() == 0:
            # Try sanitize=False as last resort
            try:
                mh = Chem.AddHs(m)
                if AllChem.EmbedMolecule(mh, useRandomCoords=True, maxAttempts=10, randomSeed=1) != 0:
                    spectators.append(Chem.MolToSmiles(m))
                    continue
            except:
                spectators.append(Chem.MolToSmiles(m))
                continue
        dn = _find_donor_simple(mh, metal)
        if dn is None:
            spectators.append(Chem.MolToSmiles(m))
            continue
        di, de, _ = dn
        neigh = [nb.GetIdx() for nb in mh.GetAtomWithIdx(di).GetNeighbors()]
        units.append((mh, di, neigh, de))

    if not units:
        return {"status": "failed", "error": "no donor atoms found"}

    # If too few donors for expected CN, add a small placeholder ligand
    expected = EXPECTED_CN.get((metal, ox), max(2, len(units)))
    while len(units) < expected and len(units) < 6:
        # Build a simple placeholder (Cl- or NH3)
        if metal in ("Au", "Ag", "Cu") and expected <= 4:
            placeholder = Chem.MolFromSmiles("[Cl-]")
            placeholder.SetProp("_is_placeholder", "1")
        else:
            placeholder = Chem.MolFromSmiles("N")
        if placeholder is None: break
        mh = _embed(placeholder, 1)
        if mh is None: break
        di, de = 0, "Cl" if "[Cl-]" in Chem.MolToSmiles(placeholder) else "N"
        neigh = [nb.GetIdx() for nb in mh.GetAtomWithIdx(di).GetNeighbors()]
        units.append((mh, di, neigh, de))
        if len(units) >= expected: break

    cn = len(units)
    if cn == 1:
        axes = _LINEAR
    elif cn == 2:
        axes = _LINEAR
    elif cn == 4:
        # Determine geometry: square planar for d8 metals, tetrahedral otherwise
        if metal in ("Au", "Ag", "Cu", "Pt", "Pd", "Ni", "Rh") and ox in (1, 2, 3):
            axes = _SQ_PLANAR
        else:
            axes = _TETRA
    elif cn == 6:
        axes = _OCTA
    else:
        # Use as many as we have
        axes = _OCTA[:min(cn, 6)] if cn <= 6 else _OCTA
    axes = [np.array(a, float) / np.linalg.norm(a) for a in axes]
    cn = min(cn, len(axes))

    # Assemble
    all_sym = [metal]
    all_pos = [np.zeros(3)]
    core_bonds = []
    ranges = []
    offset = 1
    for (mh, di, neigh, de), axis in zip(units[:cn], axes):
        conf = mh.GetConformer()
        pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])
        r = _R_DONOR.get(de, 2.1)
        moved = _align(pos, di, neigh, axis, r)
        all_sym += [a.GetSymbol() for a in mh.GetAtoms()]
        all_pos.append(moved)
        core_bonds.append((0, offset + di, r))
        ranges.append((offset, offset + mh.GetNumAtoms()))
        offset += mh.GetNumAtoms()
    all_pos = np.vstack(all_pos)

    # Set exact core-donor distances
    for i, j, d in core_bonds:
        v = all_pos[j] - all_pos[i]
        n = np.linalg.norm(v)
        if n > 1e-6: all_pos[j] = all_pos[i] + v / n * d

    # Declash
    all_pos = _declash(all_pos, all_sym, ranges, axes, passes=5)

    # Total atoms check
    if len(all_sym) > max_atoms:
        return {"status": "failed", "error": f"too many atoms: {len(all_sym)} > {max_atoms}"}

    # Charge
    charge = sum(Chem.GetFormalCharge(Chem.RemoveHs(u[0])) for u in units) + ox

    return {
        "status": "ok", "method": "simple_builder",
        "symbols": all_sym, "positions": all_pos,
        "n_atoms": len(all_sym), "total_charge": int(charge),
        "metal": metal, "cn": cn, "spectators": spectators,
    }


def write_xyz(result, path, title=""):
    s, x = result["symbols"], result["positions"]
    with open(path, "w") as f:
        f.write(f"{len(s)}\n{title} charge={result['total_charge']} "
                f"method={result['method']} cn={result['cn']}\n")
        for sym, p in zip(s, x):
            f.write(f"{sym} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    return path


# ═══════════════════════════════════════════════════════════════════
#  Main: run on all 194 failed complexes
# ═══════════════════════════════════════════════════════════════════

def main():
    import sqlite3, urllib.request
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="/root/biometaldb-3d/pilot.sqlite")
    ap.add_argument("--out", default="/root/biometaldb-3d/polynuclear/out_simple")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # Load failed IDs from VPS index.json (once)
    url = "https://mol.biometal.xyz/viewer/full/index.json"
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            data = json.load(r)
    except Exception as e:
        print(f"Cannot fetch VPS index.json: {e}")
        return
    failed = [r["id"] for r in data["records"] if r.get("status") == "no_structure"]
    print(f"Failed IDs from VPS: {len(failed)}", flush=True)

    con = sqlite3.connect(args.db)
    ok = 0
    by_reason = {}
    for cid in failed:
        row = con.execute(
            "SELECT metal, oxidation_state, smiles_ligands FROM complexes WHERE id=?",
            (cid,)).fetchone()
        if not row: continue
        metal, ox, smi = row
        try:
            r = build(metal, ox, smi, max_atoms=200)
        except Exception as e:
            r = {"status": "failed", "error": f"EXC {type(e).__name__}: {str(e)[:80]}"}
        if r["status"] == "ok":
            ok += 1
            write_xyz(r, os.path.join(args.out, f"simple_{cid}.xyz"),
                      title=f"BiometalDB #{cid} {metal}({ox})")
            if ok % 10 == 0 or ok < 5:
                print(f"  OK #{cid} {metal}({ox}) cn={r['cn']} atoms={r['n_atoms']} ({ok}/{len(failed)})", flush=True)
        else:
            by_reason.setdefault(r["error"][:60], 0)
            by_reason[r["error"][:60]] += 1
    con.close()
    print(f"\nDONE: built {ok}/{len(failed)}")
    for err, n in sorted(by_reason.items(), key=lambda x: -x[1])[:10]:
        print(f"  {n:3d}  {err}")


if __name__ == "__main__":
    main()
