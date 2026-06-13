#!/bin/bash
# Rebuild the recon3d compute env on hive (idempotent, detached-safe).
set -e
cd /root/biometaldb-3d
rm -rf arch-env mamba/pkgs 2>/dev/null || true
echo "[setup] $(date) fetching micromamba"
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=/root/biometaldb-3d/mamba
echo "[setup] creating conda env"
./bin/micromamba create -y -p /root/biometaldb-3d/arch-env -c conda-forge \
  "python=3.11" xtb-python openbabel rdkit pynauty mendeleev py3dmol numpy scipy pandas numba matplotlib
RUN="./bin/micromamba run -p /root/biometaldb-3d/arch-env"
echo "[setup] pip: architector"
$RUN pip install -q architector==0.0.10 --no-deps
echo "[setup] pip: torch cpu"
$RUN pip install -q torch --index-url https://download.pytorch.org/whl/cpu
[ -d pydentate ] || git clone --depth 1 https://github.com/hjkgrp/pydentate.git
echo "[setup] pip: pydentate + trexio + ase"
$RUN pip install -q -e ./pydentate trexio ase
echo "[setup] verifying"
$RUN python -c "
import warnings; warnings.filterwarnings('ignore')
import architector, torch, trexio
from xtb.ase.calculator import XTB
import sys; sys.path.insert(0,'recon3d')
from pydentate.pydentate_lite import pydentate_lite
import coord_oracle, assign, build, validate, trexio_writer, run_ir100
print('HIVE_ENV_READY', torch.__version__)
"
touch /root/biometaldb-3d/SETUP_DONE
echo "[setup] $(date) DONE"
