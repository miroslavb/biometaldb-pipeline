"""Deep diagnostic for #7874 with bent_template_v2 — trace assembly + xTB."""
import warnings; warnings.filterwarnings("ignore")
import sys, os, subprocess, tempfile, time, shutil
sys.path.insert(0, "/root/biometaldb-3d/polynuclear")
import numpy as np
from rdkit import Chem
import bent_template_v2 as BT

import sqlite3
con = sqlite3.connect("file:/root/biometaldb-3d/pilot.sqlite?mode=ro", uri=True)
metal, ox, smi = con.execute(
    "SELECT metal,oxidation_state,smiles_ligands FROM complexes WHERE id=7874").fetchone()
con.close()

frags = [f.strip() for f in smi.split(".") if f.strip()]
print(f"#7874: {metal}({ox})  fragments={len(frags)}")

parsed = []
for f in frags:
    m = Chem.MolFromSmiles(f) or Chem.MolFromSmiles(f, sanitize=False)
    parsed.append(m)

units = []
for i, m in enumerate(parsed):
    mh = BT._embed(m, 1)
    if mh is None:
        print(f"  frag{i} EMBED FAILED")
        continue
    dn = BT._find_core_donor(mh)
    is_met, msym, midx = BT._is_metallocene(mh)
    print(f"  frag{i}: atoms={mh.GetNumAtoms()} donor={dn[1] if dn else 'NONE'} metallocene={is_met}")
    if dn is None:
        continue
    di, de, _ = dn
    neigh = [nb.GetIdx() for nb in mh.GetAtomWithIdx(di).GetNeighbors()]
    units.append((mh, di, neigh, de))

print(f"\nUnits: {len(units)}")

prepped = []
for mh, di, neigh, de in units:
    syms, pos = BT._prep_fragment(mh, 1)
    prepped.append((syms, pos, di, neigh, de, mh))
    # Check for zero distances
    from itertools import combinations
    min_d = 999
    for a, b in combinations(range(len(syms)), 2):
        d = np.linalg.norm(pos[a] - pos[b])
        if d < 0.01:
            print(f"  ZERO-DIST in prepped: {syms[a]}{a}–{syms[b]}{b} = {d:.4f}")
        if d < min_d:
            min_d = d
    print(f"  prepped: {len(syms)} atoms min_d={min_d:.3f}Å")

# Assemble
axes = [np.array(a, float) / np.linalg.norm(a) for a in BT._GEOM_AXES[2]]
all_sym = [metal]
all_pos = [np.zeros(3)]
core_bonds = []
internal_bonds = []
ranges = []
frag_axes = []
offset = 1
tot_charge = 0

for (syms, pos, di, neigh, de, mh), axis in zip(prepped, axes):
    r = BT._R_DONOR.get(de, 2.1)
    moved = BT._align(pos, di, neigh, axis, r)
    tot_charge += Chem.GetFormalCharge(Chem.RemoveHs(mh))
    for mi, nj, d in BT._internal_metal_bonds(mh):
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

print(f"\nPre-declash:")
for a in range(len(ranges)):
    for b in range(a + 1, len(ranges)):
        d = BT._min_cross(all_pos, all_sym, range(*ranges[a]), range(*ranges[b]))
        print(f"  frag{a}-frag{b}: {d:.2f} Å")

all_pos = BT._declash(all_pos, all_sym, ranges, frag_axes, passes=5)
print(f"Post-declash:")
for a in range(len(ranges)):
    for b in range(a + 1, len(ranges)):
        d = BT._min_cross(all_pos, all_sym, range(*ranges[a]), range(*ranges[b]))
        print(f"  frag{a}-frag{b}: {d:.2f} Å")

# Check for zero distances in assembled
min_d = 999
for a, b in combinations(range(len(all_sym)), 2):
    d = np.linalg.norm(all_pos[a] - all_pos[b])
    if d < 0.01:
        print(f"  ZERO-DIST assembled: {all_sym[a]}{a}–{all_sym[b]}{b} = {d:.4f}")
    if d < min_d:
        min_d = d
print(f"Assembled min_d={min_d:.3f}Å")

with open("/tmp/diag5_assembled.xyz", "w") as f:
    f.write(f"{len(all_sym)}\n#7874 assembled v2\n")
    for s, p in zip(all_sym, all_pos):
        f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

charge = tot_charge + ox
print(f"\nCharge={charge} atoms={len(all_sym)}")

# xTB single-point
d = tempfile.mkdtemp(prefix="diag5_")
with open(os.path.join(d, "in.xyz"), "w") as f:
    f.write(f"{len(all_sym)}\n\n")
    for s, p in zip(all_sym, all_pos):
        f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

for etemp in [1500, 3000, 5000, 8000]:
    cmd = ["xtb", "in.xyz", "--gfn", "2", "--chrg", str(int(charge)),
           "--uhf", "0", "--etemp", str(etemp), "--acc", "2.0"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=120,
                       env=dict(os.environ, OMP_NUM_THREADS="6"))
    wall = time.time() - t0
    scf_ok = "SCF converged" in r.stdout or "convergence criteria satisfied" in r.stdout
    energy = None
    for line in r.stdout.splitlines():
        if "TOTAL ENERGY" in line:
            try:
                energy = float(line.split()[-2])
            except:
                pass
    print(f"  GFN2 etemp={etemp}: rc={r.returncode} wall={wall:.0f}s "
          f"SCF={'OK' if scf_ok else 'FAIL'} E={energy}")

shutil.rmtree(d, ignore_errors=True)
print("DIAG5_DONE")
