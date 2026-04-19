# ORCA on GCE — Procedure

## Instance: orca-ir-opt (34.23.233.62)
- 4 vCPU, 16GB RAM, us-east1-b
- ORCA 6.1.1: `/root/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`
- System OpenMPI 4.1.6

## MPI Parallelism (SOLVED 2026-04-18)

ORCA internally calls `mpirun -np N` WITHOUT `--allow-run-as-root` and WITHOUT `--oversubscribe`. The fix is environment variables, not a different OpenMPI version.

**Use system OpenMPI 4.1.6** — do NOT compile 4.1.8 (ABI incompatibility causes `MPI_Type_match_size` error with RIJCOSX).

```bash
ORCA_DIR=/root/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export LD_LIBRARY_PATH=$ORCA_DIR/lib:$ORCA_DIR:$LD_LIBRARY_PATH
export ORCA_TMPDIR=/dev/shm/orca_tmp
export OMPI_MCA_orte_tmpdir_base=/dev/shm/orte_sessions
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_rmaps_base_oversubscribe=1
mkdir -p $ORCA_TMPDIR /dev/shm/orte_sessions
$ORCA_DIR/orca input.inp > output.out 2>&1 &
```

## Input Rules
- `%scf ... end` — MUST close with `end`
- `%pal nprocs N end` — add BEFORE option blocks
- `* xyz charge multiplicity` — odd electrons -> multiplicity 2 (doublet)

## Performance
- 66 atoms PBE0/def2-SVP: ~45s/SCF iteration on 4 cores (vs 160s single-core)
- ~3.5x speedup
- Geometry opt: ~1-2 min/step with 4 cores, LooseOpt typically needs 25-40 steps

## Geometry Extraction from Running Job (2026-04-19)

**Trajectory file:** ORCA writes to `$ORCA_TMPDIR/<name>_trj.xyz` (RAM disk). Each geometry step is appended as a 66+2-line block.

**Extract last geometry:**
```bash
# Find line number of last "66" header in trajectory
grep -n '^66$' /dev/shm/ir_opt_loose/Ir_opt_loose_trj.xyz | tail -1
# Extract that block (68 lines: count + title + 66 atoms)
sed -n '1769,1836p' /dev/shm/ir_opt_loose/Ir_opt_loose_trj.xyz > final.xyz
```

**Even unfinished jobs produce usable geometry** — 20-28 steps of DFT optimization gives correct bond lengths and topology. Formal convergence (RMS gradient < 0.0005) is cosmetic.

## Converting XYZ → SDF/PDB with RDKit + Viewer

**Problem:** `Chem.MolFromXYZFile()` fails on ORCA output (no bond info). Must build mol manually.

**Solution:**
1. Parse XYZ manually (symbol + coords)
2. Create `RWMol`, add atoms, add conformer
3. Add bonds by distance with **different cutoffs** for metal vs organic:
   - Metal (Ir, Z=77): `(r_cov_Ir + r_cov_X) * 1.15` — stricter, avoids false Ir-C bonds
   - Organic: `(r_cov_C + r_cov_H) * 1.3`
4. Write SDF via `SDWriter`, PDB via `PDBWriter`

**Covalent radii (Å):** H=0.31, C=0.76, N=0.71, Ir=1.41

**Critical:** RDKit SDF V2000 supports max 999 atoms and 999 bonds — fine for our complexes.

## 3Dmol.js Viewer Integration

**Problem:** 3Dmol.js auto-bond detection from XYZ fails for metal complexes (no metal-specific radii). Fetch from server may fail through reverse proxy (502).

**Solution:** Embed MOL block (V2000) **directly in HTML** as JavaScript template literal:
```javascript
let moldata=`<MOL block from RDKit MolToMolBlock()>`;
v.addModel(moldata,"sdf");
```

**Style settings:**
- Ir: `{sphere:{scale:0.5,color:'#FFD700'}}` — hex color, 'gold' string not recognized
- Ligands: `stick:{radius:0.1,colorscheme:'Jmol'}` or ball+stick with `sphere:{scale:0.25}`
- Files: `/root/.hermes-agent2/biometaldb/viewer/Ir_opt_loose.html`

## Speed Optimization Suggestions

### 1. Start with LooseOpt, refine with NormalOpt
- `LooseOpt` converges in ~25-30 steps (fast, good enough for structures)
- Only use `TightOpt` for final refinement or frequency calculations
- Two-stage: LooseOpt → restart from last geometry with NormalOpt (fewer steps needed)

### 2. Use RIJCOSX + def2-SVP (already doing this)
- RIJCOSX is essential for speed with hybrid functionals
- def2-SVP is the minimum acceptable basis for geometry (not single-point energies)

### 3. Pre-optimize with xTB (GFN2-XTB) before ORCA
- xTB takes ~5-10 seconds per complex (vs ~30 min for ORCA)
- Use `xtb input.xyz --opt --gfn 2 --chrg <charge>`
- Then feed xTB-optimized geometry to ORCA as starting point
- **Caveat:** xTB GFN2 is NOT parameterized for Ir (5d metals) — bond lengths may be wrong, but topology/connectivity will be correct
- For Co, Ru, Os (3d/4d): xTB is well-parameterized and excellent as pre-optimizer

### 4. Reduce MaxIter in %scf
- `MaxIter 100` usually sufficient (vs default 300)
- Saves time per geometry step

### 5. Use /dev/shm for all I/O
- Already implemented — eliminates file-locking bottlenecks
- Critical for RIJCOSX grid calculations

### 6. Consider GFN1-FF for pre-screening
- GFN1-FF (force field) takes <1 second
- Good for checking connectivity before expensive DFT

### 7. Batch multiple complexes sequentially
- Write a script that loops over complexes, reusing the same ORCA instance
- Avoids process startup overhead (~5-10s per job)
