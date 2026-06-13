# Checkpoint — BiometalDB 3D reconstruction

**Stamp:** 2026-06-13 22:21 UTC (Sat) · branch `feat/3d-reconstruction` @ **`17e39dd`** (pushed origin)
**Live site:** https://mol.biometal.xyz/  → "🧪 3D Structure Reconstruction" + "🧬 Ligand Library"
**Slug for return:** `biometaldb-3d-reconstruction` (gbrain page `projects/biometaldb-3d-reconstruction`)

---

## What this project is
Reconstruct **3D structural geometry** (structural-formula accuracy, QM-optimization-ready) for the
**9,414 antiproliferative metal complexes** at https://mol.biometal.xyz/complexes. The DB only retained
metal-center + ligand SMILES; the old "3D MOL viewer" structures were actually RDKit 2D with a disconnected
metal → reframed as a true 3D build pipeline.

**Pipeline** = pydentate coordination oracle → priors (geometry/spin/ox) → Architector build (xtb/UFF) →
stereoisomer enumeration (diastereomers via stereo-signature dedup + Λ/Δ enantiomers via mirroring) →
validation gates → TREXIO record (h5 + TEXT) + per-compound .zip archive.

## Key locations
| What | Where |
|---|---|
| Scratch/working tree (NOT a git repo) | `/root/biometaldb-3d/` |
| micromamba env | `/root/biometaldb-3d/arch-env` (root prefix `/root/biometaldb-3d/mamba`) |
| Git repo | `/root/.hermes-agent2/biometaldb/.git` |
| **Worktree (commit here)** | `/root/biometaldb-3d/wt` on branch `feat/3d-reconstruction` |
| recon3d code in repo | `wt/scripts/recon3d/` (mirrors scratch `/root/biometaldb-3d/recon3d/`) |
| Deployed viewer root | `/root/.hermes-agent2/biometaldb/viewer/` (Caddy: mol.biometal.xyz → Flask :8502) |
| Homepage route | `/root/.hermes-agent2/biometaldb/mol_server.py` (`/` route) |
| Oracle cache | `/tmp/oracle_full.json` (pydentate predictions, 5,090 ligands) |
| Ir-100 output | `/root/biometaldb-3d/out/ir100` → deployed `/viewer/ir100/review.html` |
| Ligand library output | `/root/biometaldb-3d/out/ligands` → deployed `/viewer/ligands/index.html` |

## Hosts
- **This box** (8-core) = orchestration + deploy + Caddy/Flask. Memory-pressured; do NOT run heavy pandas/builds here.
- **hive T06** = `100.89.172.93` (56c/125G). **hive T07** = `100.121.152.8` (current full-DB runner).
- Biometal data is **NOT PII** → OK to compute on hive.
- TG notify: `/root/biometaldb-3d/hive_notify.sh` (chat_id 54300857) — **contains the bot token, MUST NOT be committed**;
  token sourced from `/root/.claude-tg-bridge/.env`.
- Detached jobs on hive: use **`systemd-run`** (setsid/nohup get killed by hive session-scope cleanup on ssh disconnect; no tmux).

---

## STATUS

### ✅ Ir-100 (done)
99/100 Ir complexes reconstructed (#4897 = `[Ir(ppy)2(dppp)]+` unbuildable by Architector AND RDKit-embed+xtb → flagged
for molSimplify; #4896 recovered at n_conf=4). **200 structures** total (Λ/Δ enantiomers for all 98 chiral). All valid;
TREXIO h5+txt + per-compound archives; review UI live with auto-load 3D (IntersectionObserver) + "✓ validate all passing" button.

### ✅ Ligand library (done, this session's main work)
**5,090 unique ligands** classified into 25 chemical groups. This session added **6 charge-resolved κ² bidentate-chelate
groups** by user request (each: edit `classify()` → `--reclassify` reuse PNGs → redeploy json → commit):

```
  175  κ²-C,N cyclometalated (anionic, −1)      ← ppy⁻-type (commit 50077af)
   61  Cyclometalated C^N (other)
    7  κ²-C,C bis-carbene chelate (neutral)     ← bis-NHC (commit 17e39dd)
   86  κ²-C,N carbene chelate (neutral)         ← pyridyl/picolyl-NHC (commit 17e39dd)
  465  N-heterocyclic carbene (NHC)             ← monodentate/other
 1243  κ²-N,N chelate (neutral diimine/polypyridyl)  ← bipy/phen (commit 7a1ef82)
  124  κ²-N,N chelate (anionic, −1)             ← amidinate/nacnac/pyridyl-pyrrolide (commit 7a1ef82)
  106  N,N-chelate (other)
```
Classification logic: gate on `pred["cn"]` (denticity from oracle) + donor-atom multiset (`syms`) + `Chem.GetFormalCharge(mol)`.
Carbene C identified as coord-C with 0 H and ≥2 N neighbours. Full group_counts dumped in `out/ligands/ligands.json`.

### 🟢 Full-DB reconstruction (IN PROGRESS on hive T07)
- **`bm-full.service`** (systemd, active/running), started Sat 2026-06-13 21:43 +03.
- **Progress: 1907 / 9414 (20.2%)** as of 00:46 +03 / 21:46 UTC.
- 1821 compounds with valid 3D, 1977 valid isomers on disk; 75 `no_structure` (flagged, like #4897); **0 hard failures**.
- Rate ~10.4 compounds/min (20 workers). **ETA ≈ 12 h → ~Sun 14 Jun ~12:45 local.** Tail is slower (big ligands: 70–136 s/cmpd vs early 19 s).
- Output: `/root/biometaldb-3d/out/full/` — `records.jsonl` (resumable), `struct/`, `img/`. Currently processing metal `Ru` (DB walked in id order).
- **On completion** `run_full_hive.sh` auto-runs post-passes (`add_enantiomers.py`, `add_trexio_text.py`, `make_archives.py`),
  touches `out/full/DONE`, and **TG-pings**.

---

## NEXT STEPS (when full run finishes)
1. **Pull** `out/full` from T07 → `tar cf - -C /root/biometaldb-3d out/full | ssh …` (the `-C` matters; a prior pull bug wiped local by running tar in wrong cwd).
2. Build a **paginated all-metal review UI** (ir100's review_ui.py is per-100; needs pagination/lazy fetch for ~9.4k × WebGL ctx limit ~16).
3. **Re-point homepage** "🧪 3D Structure Reconstruction" link from `/viewer/ir100/review.html` → the full review.
4. Re-run `ligands_library.py` (full, not `--reclassify`) if the ligand set changed — though library is DB-wide already so likely stable.
5. Update gbrain `projects/biometaldb-3d-reconstruction` + timeline.

## RESUME COMMANDS
```bash
# check full-DB progress
ssh root@100.121.152.8 'systemctl is-active bm-full.service; wc -l < /root/biometaldb-3d/out/full/records.jsonl; tail -2 /root/biometaldb-3d/out/full.log'

# reclassify ligands (reuse existing 5,090 PNGs — never re-render) + redeploy + commit
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
./bin/micromamba run -p arch-env python recon3d/ligands_library.py --db pilot.sqlite --oracle /tmp/oracle_full.json --out out/ligands --reclassify
cp out/ligands/ligands.json /root/.hermes-agent2/biometaldb/viewer/ligands/ligands.json
cd /root/biometaldb-3d/wt && cp /root/biometaldb-3d/recon3d/ligands_library.py scripts/recon3d/ && git add -A && git commit && git push origin feat/3d-reconstruction
```

## GOTCHAS
- Per-worker thread pinning required: `OMP_NUM_THREADS=1 OPENBLAS=1 MKL=1 NUMBA=1` (else 16 workers × ~30 threads → futex stall).
- Architector needs **explicit `metal_spin`** (else KeyError from mendeleev); low-spin parity = d_count%2.
- `/viewer/<dir>/` with trailing slash (no file) → 500; always link explicit `index.html`.
- Heavy embedded base64 in review HTML → broken link; keep UIs **data-driven** (fetch img/struct/json on demand).
- `hive_notify.sh` has the TG token → keep out of git.
- ssh commands with in-command `sleep` time out → keep commands short; launch `systemd-run` as the sole command.
