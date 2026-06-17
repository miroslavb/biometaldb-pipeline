# BiometalDB Pipeline

Enrichment, scoring, and web dashboard for coordination compound databases. Processes metal-organic complexes (Ru, Ir, Rh, Os, Re, Au) with cytotoxicity data from published literature. Includes a full 3D structure reconstruction pipeline producing genuine metal-coordinated geometries (9,414/9,414 — 100% coverage).

## Data Source

Database originates from **MetalCytoToxDB** (Krasnov et al., ChemRxiv 2025):

- **Zenodo**: [10.5281/zenodo.15853577](https://doi.org/10.5281/zenodo.15853577) — ML-assisted cytotoxicity database for transition metal complexes
- **Web mirror**: [biometaldb.streamlit.app](https://biometaldb.streamlit.app/) — original Streamlit interface
- **Local copy**: `data/biometaldb.sqlite` — 9,414 complexes, 26,801 measurements, 1,921 papers

## Architecture

```
biometaldb-pipeline/
├── mol_server.py              # Flask dashboard (port 8502) — main entry point
├── complexes_routes.py         # /complexes — paginated list, filters, cards, 3D viewer
├── dmpnn_routes.py             # /dmpnn — D-MPNN scoring results, version comparison
├── rdmetallics_scraper.py      # Async scraper for coordinate.rdmetallics.net API
├── app.py                      # Legacy Streamlit app (replaced by Flask)
├── chem_pipeline_lib/          # Drawing library (migrated from chem-pipeline)
│   └── coordination_draw.py    # PNG/SDF/CDXML export, PubChem fetch, square planar builder
├── scripts/
│   ├── enrich_donor_atoms.py   # RDKit donor identification (per-ligand splitting)
│   ├── enrich_rdmetallics.py   # RDMetallics CSD template matching
│   ├── enrich_papers.py        # CrossRef paper metadata enrichment
│   ├── generate_structures.py  # MOL V3000 + TUCAN identifier generation
│   ├── import_csv.py           # CSV → SQLite import
│   ├── xgboost_scorer.py       # Activity prediction (XGBoost)
│   ├── xgboost_scorer_v2.py    # Activity prediction (improved)
│   ├── batch_enrich_do.py      # Batch enrichment on DigitalOcean droplets
│   ├── dmpnn_batch_do.py       # Batch D-MPNN scoring on DO/RunPod
│   ├── setup_do.sh             # DO droplet setup script
│   └── recon3d/                # 3D reconstruction pipeline (feat/3d-reconstruction)
│       ├── run.py              # Runner: sample/ids/metal → qc.json + qc.md
│       ├── priors.py           # (metal, ox) → CN, polyhedron, multiplicity
│       ├── assign.py           # Fragment split → donor selection (CN-driven)
│       ├── build.py            # Architector + xTB retry ladder (UFF→GFN2)
│       ├── validate.py         # Gates: coordination, clash, bond lengths
│       ├── simple_builder.py   # RDKit-only fallback (no xTB), 191/194 recovered
│       ├── build_poly.py       # Polynuclear builder (ferrocene/titanocene tier)
│       ├── bent_template_v2.py # Dative→single bond fix, titanocene yield 3/9→8/9
│       ├── review_full_ui.py   # Full-DB review UI generator (index.json + review.html)
│       ├── build_trex.py       # Batch T-REX notation generator
│       ├── trexio_writer.py    # TREXIO text record writer
│       ├── add_enantiomers.py  # Λ/Δ enantiomer generation
│       ├── ligands_library.py  # Ligand template library
│       ├── coord_oracle.py     # pydentate oracle for donor assignment
│       └── polynuclear/        # Polynuclear-specific builders
│           ├── build_poly_v2.py # Metallocene split, dummy removal, P⁺, PEG truncation
│           └── bent_template.py # Bent metallocene template builder
├── patches/
│   ├── selfies_metal.py        # SELFIES library patch for metal support
│   └── rdkit_selfies_metal.py  # RDKit-SELFIES integration patch
├── viewer/                     # 3Dmol.js interactive 3D viewer
├── datasette/                  # Datasette config (CSS, metadata)
├── static/                     # Favicon
├── tests/                      # Tests (coordination_draw)
└── deploy.sh                   # Git pull + server restart
```

## Database Schema

```sql
-- Core tables
complexes       -- id, metal, oxidation_state, charge_complex, smiles_ligands,
                  -- donor_atoms, tucan, mol_quality, smiles_metal, has_mol3
measurements    -- complex_id, cell_line, ic50_dark, ic50_light, doi, year, abbreviation
papers          -- doi, title, authors, year, journal

-- D-MPNN scoring results
dmpnn_summary   -- complex_id (UNIQUE with scoring_version), metal, oxidation_state,
                  -- donor_atoms_gt, top1_donor_pred, match_score, n_candidates, scoring_version
dmpnn_results   -- complex_id, idx, ligand_smi, metal_smi, candidate_smi,
                  -- metal_ox, score, rank, is_top1, scoring_version
```

**Scoring versions**: `single-ligand` (deprecated, wrong for multi-ligand), `sequential` (correct multi-ligand), `rdmetallics.net` (author web API).

## Dependencies

```bash
# Core
rdkit-pypi>=2022.9
flask
requests
sqlite3

# Enrichment
tucan @ git+https://github.com/TUCAN-nest/TUCAN

# D-MPNN scoring
chemprop==2.2.3
rdchiral
rdmetallics @ git+https://github.com/moldagulovg/RDMetallics_coordinate.git
pandas
numpy
tqdm

# Drawing
Pillow>=9.0

# Web
datasette

# Dev
pytest
pytest-cov
```

## Installation

```bash
git clone https://github.com/miroslavb/biometaldb-pipeline.git
cd biometaldb-pipeline
pip install -e ".[dev]"

# TUCAN
pip install "tucan @ git+https://github.com/TUCAN-nest/TUCAN"

# RDMetallics + D-MPNN
pip install chemprop==2.2.3 rdchiral
pip install "git+https://github.com/moldagulovg/RDMetallics_coordinate.git"
```

## Deploy (GitHub → Production)

```bash
# 1. Edit code locally
# 2. Commit + push
git add -A && git commit -m "description" && git push origin master

# 3. On server:
cd /root/.hermes-agent2/biometaldb
bash deploy.sh
```

`deploy.sh` does: `git pull origin master` → kill old process → `nohup python3 mol_server.py` → verify.

**Never edit files directly on production.** Always use the GitHub workflow.

## Running the Dashboard

```bash
# Development
python3 mol_server.py

# Production (systemd or nohup)
nohup python3 mol_server.py > /tmp/mol_server.log 2>&1 &

# Access
# Complexes:  http://localhost:8502/complexes
# D-MPNN:     http://localhost:8502/dmpnn
# 3D Viewer:  http://localhost:8502/viewer/?id=<cid>
# API XYZ:    http://localhost:8502/api/structures/<cid>/xyz
```

## Enrichment Procedures

### 1. Donor Atom Identification (RDKit)

```bash
python3 scripts/enrich_donor_atoms.py
```

- Splits `smiles_ligands` via `Chem.GetMolFrags()` → per-ligand analysis → merged Counter
- Uses `sanitize=False` for metal-containing SMILES
- Heuristic: formal charge < 0, N/P with lone pairs, halides
- Result: JSON string `{"N": 2, "O": 1}` stored in `complexes.donor_atoms`
- **Performance**: ~1.6ms/complex, 9414 ≈ 15 seconds

### 2. RDMetallics Coordination Prediction

```bash
python3 scripts/enrich_rdmetallics.py
```

- Uses CSD SMARTS templates (13,716 patterns) to predict coordination modes
- `coordinate_ligand(ligand_mol, metal_mol)` → candidate coordination structures
- Without D-MPNN scoring, cannot select "correct" mode (returns all matches)
- First call slow (~10-30s loads pickle), subsequent ~0.5-1.5s
- Fallback for complexes with NULL donor_atoms (8/9414)

### 3. TUCAN Canonical Identifiers

```bash
python3 scripts/generate_structures.py
```

- Input: MOL V3000 files (from `data/mol3/`)
- TUCAN generates canonical strings independent of bonding concepts
- Format: `C35ClN7PRu/(1-2)(1-10)(2-33)...`
- Enables deduplication: `SELECT tucan, COUNT(*) ... GROUP BY tucan HAVING cnt > 1`
- Performance: 9414 complexes ≈ 40 seconds

### 4. Paper Metadata (CrossRef)

```bash
python3 scripts/enrich_papers.py
```

- Fetches title, authors, journal from CrossRef API by DOI
- Resume-capable: skips DOIs already in `papers` table
- Rate limit: ~50 req/s with polite User-Agent

### 5. XGBoost Activity Scoring

```bash
python3 scripts/xgboost_scorer.py
python3 scripts/xgboost_scorer_v2.py
```

- Predicts IC50 cytotoxicity from structural features
- Trained on MetalCytoToxDB measurements

## D-MPNN Scoring (RDMetallics + Chemprop)

Predicts donor atoms using pre-trained D-MPNN ensemble (5 models, Moldagulov et al., Angew. Chem. Int. Ed. 2026, e24655).

### Key Concepts

- **Source**: [moldagulovg/RDMetallics_coordinate](https://github.com/moldagulovg/RDMetallics_coordinate)
- **Models**: 15 pre-trained (3 splits × 5 ensemble), `model_FP_split` recommended
- **Multi-ligand**: sequential approach — each ligand processed separately, top-1 becomes metal for next
- **Match score**: Jaccard similarity on donor element counts (intersection/union)

### Scoring Methods

#### Method A: Local Sequential (single machine)

```bash
cd ~/RDMetallics_coordinate
source ~/.hermes-agent2/hermes-agent/venv/bin/activate
PYTHONUNBUFFERED=1 python3 batch_sequential.py --metal Ir --count 5 --workers 1
```

- ~3-4 min per complex (depends on n_ligands)
- Safe for small batches (5-10), use `--workers 1` on low-RAM machines

#### Method B: Batch on DigitalOcean Droplets

```bash
# On droplet (s-4vcpu-8gb, $0.071/hr):
python3 scripts/dmpnn_batch_do.py --offset 0 --limit 1800 --workers 4
```

- Single-CSV approach: all candidates in one CSV, 5 chemprop calls total
- ~16-28s/complex with 4 workers
- 8 droplets × 848 each = full DB in ~6.6 hours, ~$1.90

#### Method C: RunPod CPU Pod (recommended for full runs)

```bash
# Create pod (REST API, NOT GraphQL):
curl -s -X POST "https://rest.runpod.io/v1/pods" \
  -H "Authorization: Bearer $RUNPOD_TOKEN" \
  -d '{
    "name": "biometal-scorer",
    "imageName": "runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04",
    "computeType": "CPU",
    "cpuFlavorIds": ["cpu5c"],
    "vcpuCount": 32,
    "containerDiskInGb": 30,
    "ports": ["22/tcp"],
    "env": {"PUBLIC_KEY": "<biometaldb_batch.pub>"}
  }'

# Setup on pod:
scp -i ~/.ssh/biometaldb_batch -P $PORT biometaldb.sqlite root@$IP:/workspace/
scp -i ~/.ssh/biometaldb_batch -P $PORT -r rdmetallics/ root@$IP:/workspace/
scp -i ~/.ssh/biometaldb_batch -P $PORT sequential_batch_scorer.py root@$IP:/workspace/
ssh -i ~/.ssh/biometaldb_batch -p $PORT root@$IP \
  'apt-get install -y libxrender1 && pip install rdkit-pypi rdchiral'

# Run:
ssh -i ~/.ssh/biometaldb_batch -p $PORT root@$IP \
  'cd /workspace && PYTHONPATH=/workspace:$PYTHONPATH nohup python3 sequential_batch_scorer.py --workers 32 --flag sequential > scoring.log 2>&1 &'
```

- cpu5c 32-vCPU at $1.12/hr
- Full DB (9414 complexes): ~70 min, ~$1.30
- Results: 1.2% exact match, 18.1% partial, 70.1% no match, mean ~0.13

#### Method D: Author Web API (coordinate.rdmetallics.net)

```python
# rdmetallics_scraper.py — async scraper for author's deployment
python3 rdmetallics_scraper.py --workers 5 --limit 100
```

- Endpoint: `https://coordinate.rdmetallics.net/coordinate`
- SSL cert expired — use explicit `ssl_ctx` in aiohttp
- Server intermittently hangs (response times 5s to 50+s)
- ~0.1 complexes/sec with 5 workers (19+ hours for full DB)
- Scoring version: `rdmetallics.net`

#### Method E: GCE Server (Gaussian 16 / xTB for 3D optimization)

Not a scoring method per se, but used for post-scoring 3D structure refinement:

```bash
ssh -i ~/.ssh/biometaldb_batch root@34.24.168.43

# xTB optimization:
xtb input.xyz --opt --gfn 2 --chrg <charge> --uhf <uhf>

# Gaussian 16:
g16 < input.gjf
```

### Scoring Pitfalls

1. **Multi-ligand bug**: Old scripts (`batch_ru_scoring.py`, `batch_all_metals.py`) pass ALL ligands concatenated as one SMILES — WRONG. Use `sequential_batch_scorer.py`.
2. **metal_smiles column is unreliable**: May contain full complex SMILES. Reconstruct from `metal` + `oxidation_state` columns: `[{metal}+{ox}]`
3. **Use `model_FP_split`**, not `model_random_split` (authors' recommendation)
4. **NEVER overwrite main DB with scp** — download to separate file, merge with `INSERT OR REPLACE`, then cleanup
5. **`dmpnn_results` must be populated** — Flask detail page queries it; summary-only writes show "N Candidates" but empty detail page

## 3D Structure Reconstruction (recon3d)

**Status**: ✅ 9,414/9,414 — 100% coverage. Every complex in BiometalDB now has a genuine metal-coordinated 3D geometry.

### Why this exists

The previous pipeline's "3D" (`has_mol3=1`, 7,058 rows) was actually RDKit **2D** layouts with a **disconnected** metal atom (all z=0, no coordination bonds). So effectively **0 / 9,414** complexes had a genuine metal-coordinated 3D structure. This pipeline produces real ones.

### Architecture

```
recon3d/
├── priors.py           # (metal, ox) → CN, polyhedron, low-spin multiplicity
├── assign.py           # Fragment split → donor selection (CN-driven, not donor_atoms regex)
├── build.py            # Architector + xTB retry ladder (UFF → drop ligType → GFN2-xTB rescue)
├── validate.py         # Gates: metal coordinated, no clash, bond lengths in range
├── run.py              # Runner over sample/ids/metal → qc.json + qc.md
├── simple_builder.py   # RDKit-only fallback (no xTB), 191/194 recovered
├── build_poly.py       # Polynuclear builder (ferrocene/titanocene tier)
├── bent_template_v2.py # Dative→single bond fix, titanocene yield 3/9 → 8/9
├── build_trex.py       # Batch T-REX notation generator
├── trexio_writer.py    # TREXIO text record writer
├── add_enantiomers.py  # Λ/Δ enantiomer generation
├── review_full_ui.py   # Full-DB review UI (index.json + review.html)
├── ligands_library.py  # Ligand template library
├── coord_oracle.py     # pydentate oracle for donor assignment
└── polynuclear/
    ├── build_poly_v2.py # Metallocene split, dummy removal, P⁺, PEG truncation, GFN-FF
    └── bent_template.py # Bent metallocene template builder
```

### 3-Layer Fallback Strategy

The pipeline uses a progressive fallback to maximize coverage:

| Layer | Builder | Method | Coverage |
|-------|---------|--------|----------|
| **1** | `build.py` | Architector + GFN2-xTB constrained relax | 9,220/9,414 (97.9%) |
| **2** | `simple_builder.py` | RDKit 3D embedding + custom metal placement + MMFF cleanup (no xTB) | 191/194 recovered |
| **3** | `build_poly_v2.py` | Metallocene split + GFN-FF relax + spectator re-attachment | 3/3 recovered |

**Final result**: 194 → 3 → 0. All 9,414 complexes have genuine 3D structures.

### Layer 3: build_poly_v2 — Key Innovations

Handles the hardest polynuclear/sandwich-decorated complexes where Architector and simple_builder both fail:

1. **Metallocene split** — detects Co/Fe/Ru/Ni/Mn/Cr/V/Zr/Hf/Mo/W with ≥2 dative bonds to aromatic C, splits them out as spectators, leaving the real donor (Se⁻, NHC-C, P, S⁻) in the core
2. **Dummy removal** — `[n*]`, `[0*]` via `EditableMol.RemoveAtom` in REVERSE order (RDKit DistGeom crashes on dummies from `FragmentOnBonds`)
3. **P⁺ phosphonium** — added as donor atom (s=88) for cationic P ligands
4. **PEG truncation** — regex `C{4,12}O` → `CCO` for PEG-bloated SMILES (e.g., 511→409 atoms for #2746)
5. **GFN-FF priority** for polynuclears — GFN2 SCF is unstable on bimetallics; GFN-FF converges reliably
6. **Metallocene reconstruction** — builds fresh Co/Fe(cp)(cp) sandwich via RDKit instead of parsing SMILES

**Recovered complexes**:

| ID | Complex | Atoms | Method | Time |
|----|---------|-------|--------|------|
| 7118 | Au(I)–Se–cobaltocene | 43 (Au+2Se+30C+8H+2Co) | GFN-FF | 3.0s |
| 7893 | Au(I)–NHC–ferrocene | 101 (Au+56C+4N+38H+2Fe) | GFN-FF | 1.0s |
| 2746 | Ru(II)–3[P⁺]–PEG–bipy | 409 (Ru+186C+6P+6O+6N+204H) | GFN-FF | 19.2s |

### Validated Coordination Classes

- Ru/Os(II), Ir/Rh(III), Re(I) octahedral polypyridyl — Ru–N 2.04 Å
- Au(I) linear (R₃P–Au–Cl / X–Au–L) — Au–P 2.32, Au–Cl 2.26 Å
- Au(III) square planar (cyclometalated C^N) — Au–N 1.94, Au–C 1.99 Å
- Half-sandwich Ru/Os(arene), Ir/Rh(Cp*) piano-stool — Ir–C(Cp*) 2.19–2.26 Å (η⁵)

### Review UI

Interactive review interface at `viewer/full/review.html`:
- **Status filters**: ok (9,414), no_structure (0), crash_quarantine, timeout, exception
- **Metal filters**: Ru (4,834), Au (2,220), Ir (1,431), Rh (339), Os (337), Re (253)
- **Review actions**: accept ✓ / reject ✗ / flag ⚑ with localStorage persistence
- **3Dmol.js viewer**: auto-loading enlarged 3D viewers per complex
- **Export**: ⬇ export reviews as JSON

### Run

```bash
# Full pipeline (isolated micromamba env)
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
./bin/micromamba run -p ./arch-env python recon3d/run.py \
    --db pilot.sqlite --sample 60 --mode xtb --out recon3d/out/pilot

# Specific complexes
./bin/micromamba run -p ./arch-env python recon3d/run.py \
    --db pilot.sqlite --ids 1,4820,5108 --mode xtb

# Fast mode (UFF only, no xTB relax — for triage)
./bin/micromamba run -p ./arch-env python recon3d/run.py \
    --db pilot.sqlite --mode fast

# Regenerate review UI
python3 scripts/recon3d/review_full_ui.py
```

### Infrastructure

| Resource | Location | Purpose |
|----------|----------|---------|
| Pipeline runtime | `/root/biometaldb-3d/` | Isolated micromamba env (arch-env), pilot.sqlite |
| Hive (xTB) | hive3070t06 (54c/125GB) | Full-DB xTB relaxation |
| Output | `/root/biometaldb-3d/out/full/` | index.json (13.5MB), manifest.json, struct/*.xyz |
| Review UI | `viewer/full/review.html` | Interactive 3D review interface |
| T-REX | `/root/biometaldb-3d/out/full/trex.json` | Batch T-REX notation records |
| Branch | `feat/3d-reconstruction` | All recon3d code lives here |

### Known Gaps / Next Steps

1. **Ligand template library** — ~300–500 frequent ligands → curated coordList + denticity + hapticity. 5,090 distinct ligands; top-1000 templates → 47% of complexes fully covered.
2. **Stereochemistry** — cis/trans, fac/mer, Λ/Δ (currently one default isomer; enantiomers generated for chiral-at-metal).
3. **Confidence tiers** + human-review queue for ambiguous assignments.
4. **T-REX MAP recovery** — 302 complexes have ligand-only T-REX (no MAP section). Deep kekulize patch in `fix_trex_final.py` recovered 11 more; remaining 290+12 need RDKit-level fix in trex upstream.

### T-REX Notation

**Status**: ✅ 9,414/9,414 — 100% coverage.

| Type | Count | Description |
|------|-------|-------------|
| Full T-REX (with MAP) | 9,112 | `xyz_to_trex` — metal, ligands, coordination map |
| Ligand-only (no MAP) | 302 | `Metal{+ox} \| L=[ SMILES:... ]` — fallback for kekulize/timeout |

**Generation pipeline**:
- `build_trex_parallel.py` — ProcessPoolExecutor (48 workers), 15s timeout, 8,993/9,414 in 17 min on hive
- `diag_trex_fails.py` — classify failures: timeout (157), kekulize (151), polynuclear (12), ok-with-30s (97)
- `fix_trex_fails.py` — retry with SanitizeMol monkey-patch + 60s timeout → +108 recovered
- `fix_trex_final.py` — deep patch (Chem.Mol, SanitizeMol, RemoveHs) + ligand-only fallback → +301 recovered

**Output**: `/root/biometaldb-3d/out/full/trex.json` (1.8 MB, 9,414 entries)

## chem_pipeline_lib (Drawing Library)

Migrated from `miroslavb/chem-pipeline`. Provides:

- `draw_compound(name, smiles/molblock/pubchem_name)` → PNG, MOL, SDF, CDXML
- `build_square_planar('Pt', [('Cl',0),('Cl',0),('N',0),('N',0)], cis_pairs=[(0,1)])` → Pt complex with cis geometry
- `fetch_pubchem_mol('oxaliplatin')` → SDF + metadata from PubChem
- `render_png(mol, path, legend=..., bond_width=2)` → ACS-style 2D rendering
- `mol_to_cdxml(mol, path)` → ChemDraw XML export

## Infrastructure

| Service | Host | Port | Description |
|---|---|---|---|
| Dashboard (Flask) | 38.19.202.28 | 8502 | Complexes browser, D-MPNN results, 3D viewer |
| Datasette | 38.19.202.28 | 8501 | SQL query interface for raw data |
| Browser Archive | 38.19.202.28 | 8505 | Article archiving pipeline |
| Hive (xTB/Gaussian) | hive3070t06 | 22 | 54c/125GB — full-DB xTB relaxation, Gaussian 16 |
| RunPod | via API | 22 | CPU/GPU pods for batch scoring |
| DigitalOcean | via API | 22 | Droplets for batch enrichment |

### Hive Server (hive3070t06)

```bash
ssh hive3070t06
```

- 54 cores, 125 GB RAM
- Gaussian 16: `/root/g16/`
- xTB 6.7: `/usr/local/bin/xtb`
- Currently running: `TGGT_free_restart` (PID 344856, step 222/1512, E≈−10021.85)
- Disk: 89% (12G/110G free on /dev/sda4)

### VPS (38.19.202.28, biometal.xyz)

```bash
ssh root@38.19.202.28
```

- Flask dashboard: `mol_server.py` (port 8502, systemd auto-restart)
- Datasette: port 8501
- Browser Archive: port 8505
- Deploy: `git pull origin master` → systemd restart
- **No public proxy for 8502** — internal only
- `index.json` cache needs server restart after regeneration

## Data Safety Rules

1. **NEVER overwrite main DB with scp** — download to separate file, merge, cleanup
2. **NEVER edit production directly** — use GitHub workflow (edit → commit → push → deploy)
3. **SQLite concurrency**: only one writer at a time
4. **RDKit segfaults**: use `sanitize=False` everywhere in batch mode
5. **RunPod container disk**: not persistent — download results before destroying pods
6. **Split ligand SMILES**: always use `GetMolFrags()`, never analyze concatenated string

## License

Research project. Data from published literature (see `papers` table for DOIs).
