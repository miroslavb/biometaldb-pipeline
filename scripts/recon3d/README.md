# recon3d — 3D structure reconstruction for BiometalDB coordination complexes

Reconstructs **metal-coordinated 3D geometry** (QM-optimization-ready) for the
9,414 antiproliferative metal complexes in BiometalDB, from only: metal,
oxidation state, dot-separated ligand SMILES, and an (imperfect) aggregate
`donor_atoms` composition.

## Why this exists
The previous pipeline's "3D" (`has_mol3=1`, 7,058 rows) is actually RDKit **2D**
layouts with a **disconnected** metal atom (all z=0, no coordination bonds) —
see `data/mol3/complex_1.mol`. So effectively **0 / 9,414** complexes had a
genuine metal-coordinated 3D structure. This pipeline produces real ones.

## Approach
Engine: **Architector** (IBM) + **GFN2-xTB** relaxation, in an isolated
micromamba env (`/root/biometaldb-3d/arch-env`). Per complex:

1. `priors.py`   — (metal, ox) → coordination number, polyhedron, low-spin
   multiplicity; counterion set; donor-atom SMARTS; haptic (Cp/arene) detection.
2. `assign.py`   — split fragments → classify counterion / halide / haptic /
   ligand → detect donor atoms + chelate denticity → **select exactly CN donors**
   (CN driven by metal+ox, *not* the unreliable `donor_atoms` regex; `donor_atoms`
   is a composition cross-check only). Haptic η5-Cp/η6-arene occupy 3 facial sites.
3. `build.py`    — Architector build with a retry ladder (UFF → drop ligType →
   GFN2-xTB rescue); haptic complexes force xtb (UFF leaves the ring floating).
   Writes XYZ (with charge/mult), mol2, SDF.
4. `validate.py` — gates: metal coordinated, no atomic clash, M-donor bond
   lengths in range, charge/mult parity. (comp/CN match reported, not gating —
   haptic ring carbons inflate raw neighbour counts.)
5. `run.py`      — runner over a sample/ids/metal, emits `qc.json` + `qc.md`.

## Validated coordination classes (textbook bond lengths)
- Ru/Os(II), Ir/Rh(III), Re(I) octahedral polypyridyl — Ru–N 2.04 Å
- Au(I) linear (R3P–Au–Cl / X–Au–L) — Au–P 2.32, Au–Cl 2.26 Å
- Au(III) square planar (cyclometalated C^N) — Au–N 1.94, Au–C 1.99 Å
- **Half-sandwich** Ru/Os(arene), Ir/Rh(Cp*) piano-stool — Ir–C(Cp*) 2.19–2.26 Å (η5)

## Known gaps / next steps (priority order)
1. **Ligand template library** (~300–500 frequent ligands → curated coordList +
   denticity + hapticity). Needed because denticity by SMARTS+thresholds is
   whack-a-mole: β-diketonate (acac/curcumin, real O^O) vs fused ether-O
   (dioxophenanthroline, non-donor) are indistinguishable by path/strength.
   Frequency analysis: 5,090 distinct ligands; top-1000 templates → 47% of
   complexes fully covered, so templates handle the common/haptic core and the
   automated path handles the tail.
2. Confidence tiers + human-review queue for ambiguous assignments.
3. Stereochemistry (cis/trans, fac/mer, Λ/Δ) — currently one default isomer.
4. Scale run on hive3070T06 (full 9,414 xtb ≈ hours on 56 cores).

## Run
```bash
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
./bin/micromamba run -p ./arch-env python recon3d/run.py \
    --db pilot.sqlite --sample 60 --mode xtb --out recon3d/out/pilot
```
`--ids 1,4820,5108` for specific complexes; `--mode fast` for UFF (no relax;
clashes remain — not QM-ready, use only for triage).

Works on a **copy** of the DB (`pilot.sqlite`); the production DB and live
server (`mol_server.py`) are never touched. Branch: `feat/3d-reconstruction`.
