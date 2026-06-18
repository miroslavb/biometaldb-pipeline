# BiometalDB Pipeline

Enrichment, scoring, and web dashboard for coordination compound databases. Processes metal-organic complexes (Ru, Ir, Rh, Os, Re, Au) with cytotoxicity data from published literature.

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
│   └── setup_do.sh             # DO droplet setup script
├── scripts/recon3d/            # 3D structure reconstruction (branch: feat/3d-reconstruction)
│   ├── priors.py               # (metal, ox) → CN, polyhedron, spin, donor SMARTS, haptic detection
│   ├── assign.py               # Fragment split → classify → select exactly CN donors
│   ├── build.py                # Architector + retry ladder (UFF → xTB rescue), isomer enumeration
│   ├── validate.py             # Gates: metal coordinated, no clash, M-donor lengths, charge/spin
│   ├── run.py                  # Batch runner with QC report (JSON + Markdown)
│   ├── simple_builder.py       # Fallback: RDKit 3D + custom assembly (no xTB)
│   ├── polynuclear/            # Polynuclear / organometallic builder
│   │   ├── build_poly.py       # Fragment assembly + constrained xTB (no OpenBabel/Architector)
│   │   └── bent_template.py    # Idealized bent/parallel metallocene templates (Cp–M–Cp angles)
│   ├── ligands_library.py      # Curated ligand templates (coordList + denticity + hapticity)
│   ├── review_ui.py            # Human-review queue for ambiguous assignments
│   └── trexio_writer.py        # TREXIO HDF5 + text export for QC
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

-- 3D reconstruction (new columns/flags)
complexes_3d    -- complex_id, status, method, geometry, cn, charge, mult,
                  -- energy, n_atoms, xyz_path, mol2_path, sdf_path, h5_path, trexio_path,
                  -- isomer_lambda, isomer_delta, valid_lambda, valid_delta,
                  -- trex_string, confidence, flags, spectator_smiles
```

**Scoring versions**: `single-ligand` (deprecated, wrong for multi-ligand), `sequential` (correct multi-ligand), `rdmetallics.net` (author web API).

### Current DB Statistics (9,414 complexes)

| Metal | Count | % |
|-------|-------|---|
| Ru    | 4,834 | 51.3% |
| Au    | 2,220 | 23.6% |
| Ir    | 1,431 | 15.2% |
| Rh    | 339   | 3.6% |
| Os    | 337   | 3.6% |
| Re    | 253   | 2.7% |

| Oxidation State | Count |
|-----------------|-------|
| 2               | 4,904 |
| 3               | 2,642 |
| 1               | 1,858 |
| 4               | 6     |
| 5               | 4     |

| D-MPNN Scoring Version | Complexes |
|------------------------|-----------|
| `sequential`           | 7,058     |
| `rdmetallics.net`      | 363       |
| `single-ligand`        | 939       |

**3D Reconstruction (full run, branch `feat/3d-reconstruction`)**:
- **9,414/9,414 complexes built** (100% coverage)
- **20,202 structures** generated (2 isomers Λ/Δ per complex)
- **19,830 valid** (passed all validation gates)
- **353 recovered** via fallback/retry

Old `has_mol3=1` (7,058 rows) were **2D layouts with disconnected metal** — not genuine 3D. The new pipeline produces real metal-coordinated 3D geometries.

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

# 3D Reconstruction (recon3d) — isolated micromamba env at /root/biometaldb-3d/arch-env
architector @ git+https://github.com/ChemMatCAS/Architector
xtb-python @ conda-forge          # GFN2-xTB via ASE calculator
ase                               # Atomic Simulation Environment
openbabel                         # mol2 → SDF conversion
mendeleev                         # Element properties
h5py                              # HDF5 for TREXIO
trexio @ git+https://github.com/trex-co/trexio-python  # TREXIO I/O

# Dev
pytest
pytest-cov
```

### 3D Env Setup (one-time)

```bash
cd /root/biometaldb-3d
# micromamba already installed at ./bin/micromamba
./bin/micromamba create -p ./arch-env python=3.11 -y
./bin/micromamba run -p ./arch-env pip install \
    rdkit-pypi mendeleev pandas numpy tqdm h5py ase \
    "architector @ git+https://github.com/ChemMatCAS/Architector" \
    "trexio @ git+https://github.com/trex-co/trexio-python"
# xtb-python from conda-forge (supplies libxtb + Python bindings)
./bin/micromamba install -p ./arch-env -c conda-forge xtb-python -y
```

**Note**: The production server venv and DB are never touched. 3D runs work on a **copy** (`pilot.sqlite`).

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
- **Performance**: ~1.6ms/complex, 9,414 ≈ 15 seconds

### 2. RDMetallics Coordination Prediction

```bash
python3 scripts/enrich_rdmetallics.py
```

- Uses CSD SMARTS templates (13,716 patterns) to predict coordination modes
- `coordinate_ligand(ligand_mol, metal_mol)` → candidate coordination structures
- Without D-MPNN scoring, cannot select "correct" mode (returns all matches)
- First call slow (~10-30s loads pickle), subsequent ~0.5-1.5s
- Fallback for complexes with NULL donor_atoms

### 3. TUCAN Canonical Identifiers

```bash
python3 scripts/generate_structures.py
```

- Input: MOL V3000 files (from `data/mol3/`)
- TUCAN generates canonical strings independent of bonding concepts
- Format: `C35ClN7PRu/(1-2)(1-10)(2-33)...`
- Enables deduplication: `SELECT tucan, COUNT(*) ... GROUP BY tucan HAVING cnt > 1`
- Performance: 9,414 complexes ≈ 40 seconds

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
  -H "Authorization: Bearer ***" \
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
- Full DB (9,414 complexes): ~52 min, ~$0.97
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

### Scoring Pitfalls

1. **Multi-ligand bug**: Old scripts (`batch_ru_scoring.py`, `batch_all_metals.py`) pass ALL ligands concatenated as one SMILES — WRONG. Use `sequential_batch_scorer.py`.
2. **metal_smiles column is unreliable**: May contain full complex SMILES. Reconstruct from `metal` + `oxidation_state` columns: `[{metal}+{ox}]`
3. **Use `model_FP_split`**, not `model_random_split` (authors' recommendation)
4. **NEVER overwrite main DB with scp** — download to separate file, merge with `INSERT OR REPLACE`, then cleanup
5. **`dmpnn_results` must be populated** — Flask detail page queries it; summary-only writes show "N Candidates" but empty detail page

---

## 3D Structure Reconstruction Pipeline (NEW)

**Branch**: `feat/3d-reconstruction`  
**Engine**: **Architector** (IBM) + **GFN2-xTB** relaxation  
**Environment**: Isolated micromamba at `/root/biometaldb-3d/arch-env`  
**Input**: metal, oxidation_state, dot-separated ligand SMILES, aggregate `donor_atoms` (imperfect)  
**Output**: QM-optimization-ready 3D geometries (XYZ, MOL2, SDF, HDF5/TREXIO)

### Why This Exists

The previous pipeline's "3D" (`has_mol3=1`, 7,058 rows) was actually RDKit **2D** layouts with a **disconnected** metal atom (all z=0, no coordination bonds) — see `data/mol3/complex_1.mol`. So effectively **0/9,414** complexes had a genuine metal-coordinated 3D structure. This pipeline produces real ones.

### Mononuclear Pipeline (`scripts/recon3d/`)

| Stage | File | Description |
|-------|------|-------------|
| **Priors** | `priors.py` | (metal, ox) → preferred CN, polyhedron, low-spin unpaired e⁻; counterion SMILES set; donor-atom SMARTS with strengths; haptic (η⁵-Cp/η⁶-arene) detection |
| **Assignment** | `assign.py` | Split fragments → classify (ligand/halide/counterion/haptic) → detect candidate donors → **select exactly CN donors** (CN from metal+ox, NOT from unreliable `donor_atoms` regex; `donor_atoms` is cross-check only). Haptic ligands occupy 3 facial sites. Oracle (pydentate) overrides for heteroatom donors. |
| **Build** | `build.py` | Architector build with retry ladder: UFF (ligType) → UFF (no ligType) → GFN2-xTB rescue. For haptic complexes, **force xTB** (UFF leaves ring floating; xtb pulls to η⁵/η⁶). **Isomer enumeration**: generates up to N distinct stereoisomers (cis/trans, fac/mer, Λ/Δ) via symmetry sampling + stereo-signature deduplication, keeps lowest-energy per isomer. |
| **Validate** | `validate.py` | Gates: metal coordinated (≥1 M–donor bond), no atomic clash (vdW), M–donor bond lengths in range (covalent_radii × 1.3), charge/mult parity. CN/comp match reported but not gating (haptic ring carbons inflate raw neighbour counts). |
| **Run** | `run.py` | Batch runner over sample/ids/metal, emits `qc.json` + `qc.md`. Per-process timeout (default 600s for xtb). |

**Validated coordination classes (textbook bond lengths):**
- Ru/Os(II), Ir/Rh(III), Re(I) **octahedral** polypyridyl — Ru–N 2.04 Å
- Au(I) **linear** (R₃P–Au–Cl, X–Au–L) — Au–P 2.32, Au–Cl 2.26 Å
- Au(III) **square planar** (cyclometalated C^N) — Au–N 1.94, Au–C 1.99, Au–Cl 2.26 Å
- **Half-sandwich** Ru/Os(arene), Ir/Rh(Cp*) piano-stool — Ir–C(Cp*) 2.19–2.26 Å (η⁵)

### Polynuclear / Organometallic Pipeline (`scripts/recon3d/polynuclear/`)

For complexes whose ligand SMILES contain a **second metal** (titanocene/ferrocene sandwich co-ligands written with dative bonds `->`/`<-`). Architector/OpenBabel cannot handle these.

| File | Approach |
|------|----------|
| `build_poly.py` | Fragment assembly + constrained GFN2-xTB relax (OpenBabel-free). Parse fragments → find each fragment's core-metal donor (carbene C / phosphine P / thiolate S / acetylide C) → embed each fragment (RDKit ETKDG) → place around core metal at correct polyhedron (CN2=linear, CN4=tetra/square) aligning donor lone-pair axis → merge + **GFN2-xTB relax with ALL metal coordination distances FROZEN** (harmonic constraints, not hard freeze) so neither core sphere nor sandwich collapses → emit XYZ + MOL2 written directly from coordinates. |
| `bent_template.py` | **Idealised bent/parallel metallocene templates** replacing RDKit ETKDG for metallocene fragments. Bent Cp–M–Cp angles: Ti=130°, Zr/Hf=135°; Parallel sandwich: Fe/Ru/Os/Co/Ni=180°. Fixed M–Cp_centroid distances. Constrained xTB relaxation with **adaptive force constant** (softer for metallocenes: 0.5/0.8 vs 1.0/1.5) to let rings breathe. Fixes root cause of 6 titanocene failures: ETKDG produced random overlapping Cp rings that clashed on assembly and crashed xTB SCF. |

### Fallback Builder (`scripts/recon3d/simple_builder.py`)

For complexes that fail in Architector/xTB pipeline:
- **No xTB constrained relax** — just RDKit 3D embedding + custom metal-coordination placement + MMFF cleanup
- Filters haptic fragments for non-haptic metals (Au, Ag, Cu don't do η⁵/η⁶)
- Auto-adds placeholder ligands (`[Cl-]`, `NH3`) to reach expected CN for Au(I/III)
- Parses in `sanitize=False` for corrupted source donors
- Much faster and more robust for "hard" complexes (CN=2/4/6)

### Run Instructions

```bash
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
./bin/micromamba run -p ./arch-env python recon3d/run.py \
    --db pilot.sqlite --sample 60 --mode xtb --out recon3d/out/pilot
```

Options:
- `--ids 1,4820,5108` — specific complexes
- `--metal Au` — filter by metal
- `--mode fast` — UFF only (no relax; clashes remain — not QM-ready, triage only)
- `--timeout 600` — per-complex timeout (seconds)

Works on a **copy** of the DB (`pilot.sqlite`); production DB and live server (`mol_server.py`) are never touched.

### Output Structure

Per complex (`out/<cid>/`):
- `complex_<cid>.xyz` — QM input with charge/mult comment line
- `complex_<cid>.mol2` — Native connectivity including metal-coordination bonds (order 1)
- `complex_<cid>.sdf` — Via OpenBabel from mol2
- `complex_<cid>.h5` — **TREXIO HDF5** (wavefunction-ready: orbitals, integrals, geometry)
- `complex_<cid>.trexio.txt` — TREXIO text dump for inspection
- Two isomers: Λ (`_only`) and Δ (`_only_ent`) with identical energy (enantiomers)

### T-REX / TREXIO Integration

The pipeline generates **T-REX strings** and **TREXIO files** for every valid structure:

**T-REX string format** (human-readable coordination summary):
```
Au{+1} | L=[ SMILES:Cn1nnnc1[S-], SMILES:c1coc(P(c2ccco2)c2ccco2)c1 ] | MAP:{ (1:7, 2:5) }
Ru{+2} | L=[ SMILES:c1ccc(-c2ccccn2)nc1, SMILES:c1cnc2c(c1)c1c(c3cccnc32)OCCO1 ] | MAP:{ (1:13, 2:11), (1:3, 3:10) }
```

Fields:
- `Metal{oxidation}` — central metal with oxidation state
- `L=[SMILES...]` — coordinating ligands in binding order
- `MAP:{...}` — atom index mapping (ligand_atom:metal_coord_site)

**TREXIO HDF5** (`*.h5`) contains:
- `nucleus` — atomic numbers, coordinates
- `electron` — basis set, MO coefficients (if computed)
- `wave_function` — determinants, CI coefficients (extensible)
- `metadata` — method, energy, charge, multiplicity, complex_id

This enables direct feeding into quantum chemistry codes (Quantum Package, NECI, TREX).

### Known Gaps / Next Steps (Priority Order)

1. **Ligand template library** (~300–500 curated ligands → coordList + denticity + hapticity). Needed because denticity by SMARTS+thresholds is whack-a-mole: β-diketonate (acac/curcumin, real O^O) vs fused ether-O (dioxophenanthroline, non-donor) are indistinguishable by path/strength. Frequency analysis: 5,090 distinct ligands; top-1000 templates → 47% of complexes fully covered.
2. **Confidence tiers + human-review queue** for ambiguous assignments (don't fabricate connectivity the literature SMILES doesn't determine). Partial UI in `review_ui.py`.
3. **Stereochemistry** (cis/trans, fac/mer, Λ/Δ) — currently one default isomer per geometry; isomer enumeration exists but needs stereochemical assignment.
4. **Persist real 3D back to main DB** under new columns/flags; wire 3Dmol.js viewer to them; keep old 2D fields separate.
5. **Scale production run on hive3070T06** (full 9,414 × xtb ≈ hours on 56 cores). Current box (hermes, 8 cores, memory-pressured) is unsuitable — background jobs killed on session resets.

See `scripts/recon3d/README.md` for detailed run instructions.

---

## chem_pipeline_lib (Drawing Library)

Migrated from `miroslavb/chem-pipeline`. Provides:

- `draw_compound(name, smiles/molblock/pubchem_name)` → PNG, MOL, SDF, CDXML
- `build_square_planar('Pt', [('Cl',0),('Cl',0),('N',0),('N',0)], cis_pairs=[(0,1)])` → Pt complex with cis geometry
- `fetch_pubchem_mol('oxaliplatin')` → SDF + metadata from PubChem
- `render_png(mol, path, legend=..., bond_width=2)` → ACS-style 2D rendering
- `mol_to_cdxml(mol, path)` → ChemDraw XML export

## Infrastructure

| Service | Host | Port | Description |
|---------|------|------|-------------|
| Dashboard (Flask) | 38.19.202.28 | 8502 | Complexes browser, D-MPNN results, 3D viewer |
| Datasette | 38.19.202.28 | 8501 | SQL query interface for raw data |
| GCE (Gaussian) | 34.24.168.43 | 22 | Gaussian 16, xTB, ORCA calculations |
| RunPod | via API | 22 | CPU/GPU pods for batch scoring |
| DigitalOcean | via API | 22 | Droplets for batch enrichment |

### GCE Server (34.24.168.43, gaussian-lideprts)

```bash
ssh -i ~/.ssh/biometaldb_batch root@34.24.168.43
```

- Gaussian 16: `/root/g16/`
- xTB 6.7: `/usr/local/bin/xtb`
- ORCA 6.1.1: `/root/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`
- MetalloGen: `/root/MetalloGen`
- XYZ viewer: port 8765

## Data Safety Rules

1. **NEVER overwrite main DB with scp** — download to separate file, merge, cleanup
2. **NEVER edit production directly** — use GitHub workflow (edit → commit → push → deploy)
3. **SQLite concurrency**: only one writer at a time
4. **RDKit segfaults**: use `sanitize=False` everywhere in batch mode
5. **RunPod container disk**: not persistent — download results before destroying pods
6. **Split ligand SMILES**: always use `GetMolFrags()`, never analyze concatenated string

## License

Research project. Data from published literature (see `papers` table for DOIs).