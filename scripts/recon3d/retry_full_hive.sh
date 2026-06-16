#!/bin/bash
# Targeted re-run of the failed tail (353 complexes) on hive T07.
# Detached-safe: launch with `setsid bash retry_full_hive.sh &`. Resumable
# (retry_tail reloads out/full/retry.jsonl). Merges into out/full/manifest.json
# (preserving the enrichments of the 9061 already-built complexes), enriches the
# freshly recovered ones, and tags the unbuildable with a clear fail_reason.
cd /root/biometaldb-3d
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
RUN="./bin/micromamba run -p /root/biometaldb-3d/arch-env"

rm -f out/full/RETRY_DONE out/full/RETRY_FAIL
echo "===== retry_full_hive start $(date -u) =====" >> out/retry_tail.log
$RUN python recon3d/retry_tail.py --db pilot.sqlite --out out/full \
     --workers 24 --timeout 900 >> out/retry_tail.log 2>&1
rc=$?
echo "===== retry_tail rc=$rc $(date -u) =====" >> out/retry_tail.log
if [ $rc -eq 0 ]; then
  touch out/full/RETRY_DONE
else
  touch out/full/RETRY_FAIL
fi
