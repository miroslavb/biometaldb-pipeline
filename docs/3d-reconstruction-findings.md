# 3D Structure Reconstruction — Research Findings (2026-06-18, updated)

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

## Full Run Results (2026-06-17, branch feat/3d-reconstruction)

| Metric | Value |
|--------|-------|
| Total complexes | 9,414 |
| Built (geometry produced) | 9,414 (100%) |
| Structures generated (2 isomers Λ/Δ) | 20,202 |
| Valid (passed all gates) | 19,830 (98.2% of built) |
| Recovered via fallback/retry | 353 |
| By metal: Ru 4,834 · Au 2,220 · Ir 1,431 · Rh 339 · Os 337 · Re 253 |

All structures: XYZ, MOL2, SDF, TREXIO HDF5 (+ text dump), T-REX string.

## Approach (pipeline: scripts/recon3d/)

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

## Polynuclear / Organometallic complexes (scripts/recon3d/polynuclear/)

For complexes with a **second metal in ligand SMILES** (titanocene/ferrocene
sandwich co-ligands written with dative bonds `->`/`<-`), Architector/OpenBabel fail.

| Module | Solution |
|--------|----------|
| `build_poly.py` | Fragment assembly + constrained GFN2-xTB relax (OpenBabel-free). Parse fragments → find each fragment's donor to core metal (carbene C / phosphine P / thiolate S) → embed per fragment → place around core at correct polyhedron (CN2=linear, CN4=tetra/sq) → merge + **GFN2-xTB with ALL M-coordination distances FROZEN** (harmonic constraints; soft springs, not hard freeze) so core sphere and internal sandwich both held. |
| `bent_template.py` | **Idealised bent/parallel metallocene templates** replacing ETKDG. Bent Cp–M–Cp: Ti=130°, Zr/Hf=135°; Parallel: Fe/Ru/Os/Co/Ni=180°. Fixed M–Cp_centroid distances. **Adaptive force constant** for xTB: softer for metallocenes (0.5/0.8 vs 1.0/1.5) to let rings breathe. Fixes 6 titanocene failures caused by ETKDG's random overlapping Cp rings crashing xTB SCF. |

## T-REX / TREXIO Output (NEW)

Every valid structure emits:

**T-REX string** (human-readable coordination summary):
```
Au{+1} | L=[ SMILES:Cn1nnnc1[S-], SMILES:c1coc(P(c2ccco2)c2ccco2)c1 ] | MAP:{ (1:7, 2:5) }
Ru{+2} | L=[ SMILES:c1ccc(-c2ccccn2)nc1, SMILES:c1cnc2c(c1)c1c(c3cccnc32)OCCO1 ] | MAP:{ (1:13, 2:11), (1:3, 3:10) }
```
Fields: `Metal{ox}`, `L=[coordinating ligands in binding order]`, `MAP:{ligand_atom:metal_site}`.

**TREXIO HDF5** (`complex_<cid>.h5`) — wavefunction-ready:
- `nucleus` — Z, coordinates
- `electron` — basis set, MO coefficients (extensible)
- `wave_function` — determinants, CI coefficients (extensible)
- `metadata` — method, energy, charge, mult, complex_id

Plus text dump `complex_<cid>.trexio.txt` for inspection.

## Fallback Builder (scripts/recon3d/simple_builder.py)

For complexes failing Architector/xTB: RDKit 3D + custom metal-centered assembly + MMFF cleanup (no xTB). Filters haptics for Au/Ag/Cu, auto-adds `[Cl-]`/`NH3` placeholders to reach CN, parses `sanitize=False`. Covers the "hard" 194 previously failed complexes.

## Known gaps / next steps (priority order)
1. **Ligand template library** (~300–500 curated ligands → coordList + denticity + hapticity) — closes denticity tail + locks haptic/common cases.
2. **Confidence tiers + human-review queue** for ambiguous assignments (don't fabricate connectivity the literature SMILES doesn't determine). Partial UI in `review_ui.py`.
3. **Stereochemistry** (cis/trans, fac/mer, Λ/Δ) — currently one default isomer per geometry; isomer enumeration exists but needs stereochemical assignment.
4. **Persist real 3D back to main DB** under new columns/flags; wire 3Dmol.js viewer to them; keep old 2D fields separate.
5. **Scale run on hive3070T06** (full 9,414 × xtb ≈ hours on 56 cores). This box (hermes, 8 cores, memory-pressured) is unsuitable — background jobs killed on session resets.

See `scripts/recon3d/README.md` for run instructions. Branch: `feat/3d-reconstruction`.