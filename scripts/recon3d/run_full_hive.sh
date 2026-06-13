#!/bin/bash
# Full-database reconstruction on hive: all ~9.4k compounds + post-passes + TG notify.
# Resumable (run_ir100 reloads out/full/records.jsonl); re-run to continue.
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
RUN="./bin/micromamba run -p /root/biometaldb-3d/arch-env"
N(){ bash /root/biometaldb-3d/hive_notify.sh "$1" >/dev/null 2>&1; }

N "▶️ biometaldb FULL run starting on $(hostname) — all ~9.4k compounds (resumable)"
if $RUN python recon3d/run_ir100.py --db pilot.sqlite --all --workers 20 --out out/full >> out/full.log 2>&1; then
  N "… build done, post-processing (enantiomers / TREXIO-text / archives)"
  $RUN python recon3d/add_enantiomers.py out/full >> out/full.log 2>&1
  $RUN python recon3d/add_trexio_text.py out/full >> out/full.log 2>&1
  $RUN python recon3d/make_archives.py   out/full >> out/full.log 2>&1
  touch out/full/DONE
  N "✅ biometaldb FULL run DONE — $(grep -E '^\[run\] done' out/full.log | tail -1)"
else
  N "❌ biometaldb FULL run FAILED — $(tail -1 out/full.log)"
fi
