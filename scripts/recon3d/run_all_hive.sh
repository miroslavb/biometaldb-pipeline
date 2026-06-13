#!/bin/bash
# Full detached pipeline on hive: build env -> run Ir-100 -> Telegram notify.
cd /root/biometaldb-3d
rm -f SETUP_DONE RUN_DONE setup.log run.log
N() { bash /root/biometaldb-3d/hive_notify.sh "$1" >/dev/null 2>&1; }

if ! bash setup_hive.sh > setup.log 2>&1; then
  N "❌ biometaldb: hive env setup FAILED — $(tail -1 setup.log)"
  exit 1
fi
N "✅ biometaldb: hive env ready — starting Ir-100 reconstruction (oracle+xtb+isomers)"

export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
RUN="./bin/micromamba run -p /root/biometaldb-3d/arch-env"
if $RUN python recon3d/run_ir100.py --db pilot.sqlite --limit 100 --workers 16 \
      --out out/ir100 > run.log 2>&1; then
  touch RUN_DONE
  N "✅ biometaldb Ir-100 DONE on hive — $(tail -1 run.log)"
else
  N "❌ biometaldb Ir-100 run FAILED — $(tail -1 run.log)"
fi
