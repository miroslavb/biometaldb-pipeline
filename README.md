# BiometalDB Pipeline

Enrichment, scoring, and web dashboard for coordination compound databases. Processes metal-organic complexes (Ru, Ir, Rh, Os, Re) with cytotoxicity data from published literature.

## Data Source

Database originates from **MetalCytoToxDB** (Krasnov et al., ChemRxiv 2025):

- **Zenodo**: [10.5281/zenodo.15853577](https://doi.org/10.5281/zenodo.15853577) — ML-assisted cytotoxicity database for transition metal complexes
- **Web mirror**: [biometaldb.streamlit.app](https://biometaldb.streamlit.app/) — original Streamlit interface
- **Local copy**: `data/biometaldb.sqlite` — 7,058 complexes, 26,801 measurements, 1,921 papers

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
- **Performance**: ~1.6ms/complex, 7058 ≈ 11 seconds

### 2. RDMetallics Coordination Prediction

```bash
python3 scripts/enrich_rdmetallics.py
```

- Uses CSD SMARTS templates (13,716 patterns) to predict coordination modes
- `coordinate_ligand(ligand_mol, metal_mol)` → candidate coordination structures
- Without D-MPNN scoring, cannot select "correct" mode (returns all matches)
- First call slow (~10-30s loads pickle), subsequent ~0.5-1.5s
- Fallback for complexes with NULL donor_atoms (8/7058)

### 3. TUCAN Canonical Identifiers

```bash
python3 scripts/generate_structures.py
```

- Input: MOL V3000 files (from `data/mol3/`)
- TUCAN generates canonical strings independent of bonding concepts
- Format: `C35ClN7PRu/(1-2)(1-10)(2-33)...`
- Enables deduplication: `SELECT tucan, COUNT(*) ... GROUP BY tucan HAVING cnt > 1`
- Performance: 7058 complexes ≈ 30 seconds

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
- Full DB (7058 complexes): ~52 min, ~$0.97
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

## 3D Structure Generation

### MetalloGen Pipeline

Used for generating 3D conformers of metal complexes for docking.

```bash
# MetalloGen installed at: /root/MetalloGen
# ORCA: /root/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca
# xTB: /usr/local/bin/xtb (v6.7.0)
```

**Approach** (hybrid, handles MetalloGen failures):
1. Convert ligand SMILES + metal data → MetalloGen m-SMILES (e.g., `[N:1]1CCCCC1.[c-:2]1ccccc1.[Ir+3]6_octahedral`)
2. Attempt MetalloGen `TMCGenerator` 3D embedding
3. Fallback: RDKit embed ligands → geometric metal placement at donor centroid → xTB optimize
4. Output: SDF, PDB, XYZ files in `/root/ir_docking_structures/`

**Calculators**:
- **ORCA** (GFN2-XTB): accurate but single-core (MPI disabled), use `LC_ALL=C` for locale issues
- **xTB** (GFN2-xTB): faster, works directly from CLI
- **Gaussian 16**: available on GCE server (34.24.168.43), used for production QM calculations

### Generated Structures

5 Ir(III) complexes with full 3D data:

| ID | SDF | PDB | xTB opt |
|---|---|---|---|
| 4852 | ✅ | ✅ | ✅ |
| 4914 | ✅ | ✅ | ✅ |
| 5554 | ✅ | ✅ | ✅ |
| 5555 | ✅ | ✅ | ✅ |
| 5556 | ✅ | ✅ | ✅ |

Files at: `/root/ir_docking_structures/`

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
