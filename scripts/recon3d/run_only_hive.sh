#!/bin/bash
# Run the Ir-100 reconstruction (env already built) + Telegram notify. Detached-safe.
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
# pin each worker to 1 thread — 16 workers x many OMP/numba threads = lock thrash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 MKL_DYNAMIC=FALSE
RUN="./bin/micromamba run -p /root/biometaldb-3d/arch-env"
N(){ bash /root/biometaldb-3d/hive_notify.sh "$1" >/dev/null 2>&1; }
rm -f RUN_DONE run.log
N "▶️ biometaldb: starting Ir-100 reconstruction on $(hostname) (oracle+xtb+isomers)"
if $RUN python recon3d/run_ir100.py --db pilot.sqlite --limit 100 --workers 16 --out out/ir100 > run.log 2>&1; then
  touch RUN_DONE
  N "✅ biometaldb Ir-100 DONE — $(tail -1 run.log)"
else
  N "❌ biometaldb Ir-100 FAILED — $(tail -1 run.log)"
fi
