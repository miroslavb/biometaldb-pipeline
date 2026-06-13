# 3D Structure Reconstruction — Research Findings (2026-06-13)

Reconstructing QM-optimization-ready 3D geometry for the 9,414 antiproliferative
metal complexes in BiometalDB, from only: metal, oxidation state, dot-separated
ligand SMILES, and an aggregate `donor_atoms` composition.

## Key finding: the existing "3D" is not 3D
The DB flags 7,058 rows `has_mol3=1 / mol_quality='good'`, but the stored
`data/mol3/complex_*.mol` files are RDKit **2D** layouts (header says `2D`, every
z-coordinate = 0) with the metal inserted as a **disconnected** fragment
(`CHG=2 VAL=-1`, zero coordination bonds) — including the "validated" #4820.
`scripts/generate_structures.py` calls `AllChem.Compute2DCoords` on
`[Ru+2].lig1.lig2`. So `has_mol3=1` means "a file was written", not "valid 3D".
**Effectively 0 / 9,414 complexes had a genuine metal-coordinated 3D structure.**
Only a handful of hand-built `viewer/*.xyz|mol` are real 3D.

The source dataset (MetalCytoToxDB.csv) carries no coordination ground truth
either — just `SMILES_Ligands`, `Counterion`, metal/ox/charge. Donor-site
assignment, denticity, hapticity, and coordination geometry must be **inferred**.
That inference is the scientific core of the task.

## Approach (new pipeline: scripts/recon3d/)
Engine: **Architector** (IBM, github.com/ChemMatCAS/Architector) + **GFN2-xTB**
relaxation, in an isolated micromamba env at `/root/biometaldb-3d/arch-env`
(no conda available; `xtb-python` from conda-forge supplies libxtb). The
production server venv and DB are never touched.

Per complex: `priors` (metal+ox → CN / polyhedron / low-spin multiplicity;
counterion set; donor SMARTS; haptic detection) → `assign` (classify fragments;
detect donors + chelate denticity; **select exactly CN donors**, CN from metal+ox)
→ `build` (Architector + retry ladder; xtb for haptic) → `validate` (metal
coordinated, no clash, M-donor bond lengths, charge/mult) → `run` (+ QC report).

### Engineering lessons (de-risking)
- **CN must come from metal+oxidation state, not `donor_atoms`.** The DB
  `donor_atoms` was produced by a fragile regex and frequently miscounts (e.g.
  Au(I) listed with 9 donors; should be 2). Driving CN from (metal, ox) is
  reliable for this 4d/5d set; `donor_atoms` is a composition cross-check only.
- **Low-spin parity rule** for 4d/5d: unpaired = d-count % 2 (d6/d8/d10 → 0).
  Required to dodge an Architector `KeyError 'Au'` (mendeleev lazy-attr) — pass
  `metal_spin` explicitly.
- **Architector quirks:** needs explicit `coordList` (no auto-donor for arbitrary
  SMILES); infers `ligType` from `coordList` but slowly (pass mono/bi_cis to speed
  up); returns empty for bulky low-CN (e.g. (PPh3)2Au+) under UFF → rescue with
  GFN2-xTB relax (retry ladder).
- **Half-sandwich** (η6-arene, η5-Cp/Cp*): pass all ring carbons as `coordList`;
  the ligand occupies **3 facial sites**; **xtb relaxation is required** (UFF
  leaves the ring floating; xtb pulls it to η5/η6). [Cp*Ir(bpy)Cl]+ → Ir–C 2.19–2.26 Å.
- **Denticity is whack-a-mole with thresholds.** β-diketonate (acac/curcumin:
  anionic-O + carbonyl-O, real 6-membered O^O) vs fused ether-O
  (dioxophenanthroline: non-donor) are indistinguishable by bond-path length +
  donor strength. Chelate-ring window 5–6-membered (path 2–4) + strongest-pocket
  + weak-tag-along drop handles most, but the genuine fix is a **ligand template
  library** (which is why molSimplify/Architector ship curated libraries).

## Validated coordination classes (textbook bond lengths)
- Ru/Os(II), Ir/Rh(III), Re(I) **octahedral** polypyridyl — Ru–N 2.04 Å
- Au(I) **linear** (R3P–Au–Cl, X–Au–L) — Au–P 2.32, Au–Cl 2.26 Å
- Au(III) **square planar** cyclometalated C^N — Au–N 1.94, Au–C 1.99, Au–Cl 2.26 Å
- **Half-sandwich** Ru/Os(arene), Ir/Rh(Cp*) piano-stool — Ir–C(Cp*) 2.19–2.26 Å (η5)

## Ligand-frequency / template-coverage analysis
9,414 complexes → **5,090 distinct ligand fragments**. Top by # complexes:
`[Cl-]` 4180 · p-cymene 1964 · bpy 677 · PPh3 616 · Cp* 529 · CO 435 ·
cyclometalated ppy 418 · dppz 414 · benzene 238 · PTA 221 · Cp 120 · dppe 114 · NHC 103.
Templating the top-K ligands fully covers: top-100 → 8%, top-300 → 21%,
top-500 → 30%, top-1000 → 47% of complexes. **Half-sandwich (cymene+Cp*+arene+Cp)
is the #2 motif after chloride** and is now handled.

## Known gaps / next steps (priority order)
1. **Ligand template library** (~300–500 curated ligands → coordList + denticity
   + hapticity) — closes the denticity tail + locks haptic/common cases.
2. Confidence tiers + human-review queue for ambiguous assignments (don't
   fabricate connectivity the literature SMILES doesn't determine).
3. Stereochemistry (cis/trans, fac/mer, Λ/Δ) — currently one default isomer.
4. **Scale run on hive3070T06** (full 9,414 × xtb ≈ hours on 56 cores). This box
   (hermes, 8 cores, memory-pressured) is unsuitable — background jobs are killed
   on session resets.
5. Persist real 3D (MOL/SDF/XYZ + regenerated TUCAN from correct connectivity)
   back to the DB under new columns/flags; wire the 3Dmol.js viewer to them;
   keep the old 2D fields separate.

See `scripts/recon3d/README.md` for run instructions. Branch: `feat/3d-reconstruction`.
