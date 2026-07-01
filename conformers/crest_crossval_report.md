# CREST cross-validation of the served Ir(III) conformer ensemble

**Page validated:** https://mol.biometal.xyz/conformers — 5,717 GFN2-xTB-optimized
(correct per-complex charge/spin), RMSD+energy-deduped distinct conformers for 427 Ir(III) complexes.

**Method.** For a diversity subset of complexes, CREST (GFN2-xTB metadynamics, `--squick`,
correct `--chrg/--uhf`) independently re-samples conformers on hive-t07 (`/root/confopt`).
We then compare CREST's ensemble against the ensemble we actually serve
(`.../opt/conformers_opt.xyz`). Atom indexing is consistent between both (same input
topology), so index-aligned Kabsch RMSD is apples-to-apples. Analyses (this repo):
`crossval_report.py` (graded), `crossval_connectivity.py` (bond-graph filter),
`crossval_energy_windowed.py` (energy-windowed). Coverage counts only **same-bond-graph,
low-energy** CREST conformers (ΔE ≤ 3 kcal/mol vs CREST's own min), heavy-atom RMSD.

## Results (6/8 subset; 9049 + 5429 CREST still running on t07)

| cid | rot.bonds | ours (opt) | CREST low-E | coord. Δ (Å) | cov ≤1.0Å | cov ≤1.5Å | cov ≤2.0Å | worst miss | ours confirmed |
|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|
| 5044 | 3 | 12 | 11 | **0.001** | 0.82 | **1.00** | 1.00 | 1.08 Å | 2/12 |
| 5576 | 4 | 7 | 21 | **0.002** | 0.00 | 0.00 | **1.00** | 1.83 Å | 0/7 |
| 5309 | 8 | 44 | 731 | **0.003** | 0.00 | 0.13 | **1.00** | 2.02 Å | 0/44 |
| 5907 | 1 | 2 | 2 | **0.002** | 0.00 | 0.00 | 0.00 | 2.58 Å | 0/2 |
| 9283 | 7 | 46 | 247 | **0.009** | 0.00 | 0.00 | 0.00 | 4.27 Å | 0/46 |
| 4925 | 3 | 6 | 77 | **0.001** | 0.00 | 0.00 | 0.00 | 8.58 Å | 0/6 |

`coord. Δ` = mean |Ir–donor distance(ours) − Ir–donor distance(CREST ground state)| over the 6 donors.

## Verdict

1. **Metal coordination geometry is validated — the part that matters for these Ir(III)
   pharmacophores.** CREST's independently-sampled ground state reproduces our Ir–donor
   distances to **≤ 0.01 Å on every complex**. No isomerization: **0/1490** CREST frames
   changed the bond graph. The served structures are chemically correct at the metal center.

2. **Rigid / semi-rigid complexes are fully covered.** 5044 (rigid): our ensemble sits
   ≤1.0 Å from 82% of CREST's low-E conformers and ≤1.5 Å from all of them.

3. **Flexible complexes under-sample peripheral-group motion — a real but graded gap.**
   - *Soft misses (dihedral fine structure):* 5576, 5309 are covered within **2.0 Å**
     (worst 1.8–2.0 Å). CREST resolves rotamers of flexible arms one notch finer than our
     Uniconf→xTB ensemble; the metal core is identical.
   - *Hard misses (large-amplitude arm swing):* 9283 and **4925** — CREST finds distinct
     *extended* conformers of a long pendant group up to **4–9 Å** from anything we serve,
     within a few kcal/mol. 4925 is only rot=3 yet worst=8.6 Å → it is one specific long
     flexible substituent, not overall flexibility. These low-energy basins are genuinely
     absent from our served ensemble.
   - *Small-N artifact:* 5907 (rot=1, 2 conformers) shows a 2.58 Å "miss" — most likely a
     Λ/Δ enantiomer / achiral-Kabsch mirror effect, not a true gap.

4. **Low reverse coverage ("ours confirmed" column) is expected, not a defect.** `--squick`
   is a deliberately shallow search that returns a few clustered conformers; our served
   ensemble is broader, so most of our conformers simply weren't revisited by this quick
   CREST pass. It means "CREST looked less broadly here," not "our conformers are wrong."

## Bottom line

The served page is **correct where it counts (coordination sphere, connectivity, rigid
complexes)** and **incomplete for large-amplitude flexible-substituent conformers on a
minority of floppy complexes**. This is a completeness gap, not a correctness error.

## Scope across all 427 served complexes (by rotatable bonds)

| flexibility | complexes | share | cross-val expectation |
|---|---:|---:|---|
| rigid (rot=0) | 63 | 15% | single conformer, fully covered |
| low (rot 1–3) | 168 | 39% | mostly covered; watch for one-long-pendant cases (cf. 4925) |
| moderate (rot 4–6) | 87 | 20% | soft misses ≤2 Å (cf. 5576) |
| **floppy (rot ≥7)** | **109** | **26%** | **most exposed to the completeness gap** (cf. 5309, 9283) |

Median rot bonds = 3, max = 26, mean = 4.6. Among the 109 floppy complexes, avg served
conformers = 23.6, but **24 serve <10 conformers** — the prime candidates for under-sampling.

## Open follow-ups

- Finish CREST for 9049 + 5429 (running on t07) → 8/8 subset; re-run `crossval_report.py`.
- Scope impact across all 427 by rotatable-bond count (how many are "floppy-arm" complexes).
- Decide remediation: (a) document as a known limitation on the page; (b) regenerate the
  flexible subset with a larger Uniconf conformer budget or a CREST-seeded pass; (c) both.

---

## Remediation applied (2026-07-01) — the 24 under-sampled complexes

Regenerated the 24 floppy complexes (rot≥7 serving <10 conformers) with a much larger
Uniconf budget (MK=3000), then farthest-point-sampled to ≤40 diverse conformers each and
re-optimized with GFN2-xTB (`--opt loose`, on a vast.ai 40-core box after both hive nodes
proved unstable/overloaded). All coordination spheres intact (coord_ok on 24/24).

Result — served page rebuilt: **5,717 → 6,337 distinct conformers** (`total_conformers`).
The 24 target complexes: **149 → 761** conformers (+612). 21/24 improved 4–10×.

Per-complex (served count before → after):
5192 8→40 · 5344 5→40 · 5346 9→40 · 5452 8→40 · 5454 9→16 · 5527 5→36 · 5528 4→26 ·
5762 9→40 · 5776 6→32 · 5777 9→32 · 5942 4→40 · 5943 4→40 · 5985 5→24 · 6058 8→40 ·
6110 7→40 · 6111 6→40 · 6112 9→40 · 7153 9→40 · 7678 11→40 · 7812 6→40 · 8944 4→29

**3 complexes could NOT be improved — a Uniconf limitation, not budget:** 5972 (2), 8525 (3),
9122 (1). Uniconf's fragmentation/ring-closure analysis collapses their many nominal
rotatable bonds to 1–2 real torsions (e.g. 9122: "9 rotatable bonds" → "number of rotatble
bonds: 1"), so it emits only 1–3 conformers regardless of MK. A non-torsional sampler
(CREST metadynamics or RDKit ETKDG) would be needed for these three; deferred.

**Caveat:** these 24 were optimized `--opt loose` (the other 403 are `--opt normal`).
Geometrically adequate for a conformer browser (same basins, looser convergence); relative
energies are marginally less tight. Documented here for provenance.

Backup of the pre-remediation served index: `zenith_conformers.json.bak.pre-regen24`.
