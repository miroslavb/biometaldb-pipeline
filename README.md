# BiometalDB Pipeline

Enrichment and analysis pipeline for coordination compound databases (Ru, Ir, Rh, Os, Re complexes).

## Components

### Scripts

- `scripts/enrich_donor_atoms.py` — Identify donor atoms from SMILES via RDKit. Splits multi-ligand SMILES into individual fragments before analysis.
- `scripts/enrich_rdmetallics.py` — Donor atom prediction via RDMetallics (CSD template matching).
- `scripts/enrich_papers.py` — Paper metadata enrichment via CrossRef API.
- `scripts/import_csv.py` — Import complexes from CSV into SQLite.
- `scripts/generate_structures.py` — Generate MOL V3000 files and TUCAN identifiers.
- `scripts/xgboost_scorer.py` / `xgboost_scorer_v2.py` — Activity scoring models.

### Web UI

- `mol_server.py` — Flask server for 2D structure rendering and MOL file downloads (port 8502).
- `app.py` — Streamlit/Datasette web interface.
- `datasette/` — Datasette configuration (CSS, metadata).

### Patches

- `patches/` — Patches for external libraries (group-selfies metal support).

## Database

SQLite database with 7058 complexes and 26801 measurements. Schema:

- `complexes` — id, metal, oxidation_state, smiles_ligands, donor_atoms, tucan, mol_quality
- `measurements` — complex_id, cell_line, ic50_dark, ic50_light, doi
- `papers` — doi, title, authors

## Key Pitfalls

### Ligand Splitting

`smiles_ligands` contains ALL ligands as one dot-separated SMILES string. Always split via `Chem.GetMolFrags()` before per-ligand analysis. Processing the whole string gives wrong donor counts.

### RDKit + Coordination Compounds

- Use `sanitize=False` for metal-containing SMILES
- `GetTotalNumHs()` crashes on unsanitized mols — wrap in try/except
- Metal is NOT in smiles_ligands (separate column)

### TUCAN

- Input: MOL V3000 files
- API: `graph_from_molfile_text()` → `canonicalize_molecule()` → `serialize_molecule()`
- Deduplication: same TUCAN string = same structure

## Requirements

```
rdkit
tucan @ git+https://github.com/TUCAN-nest/TUCAN
rdmetallics
datasette
flask
requests
```
