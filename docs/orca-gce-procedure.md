# ORCA on GCE — Procedure

## Instance: orca-ir-opt (34.23.233.62)
- 4 vCPU, 16GB RAM, us-east1-b
- ORCA 6.1.1: `/root/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`
- System OpenMPI 4.1.6

## MPI Parallelism (SOLVED 2026-04-18)

ORCA internally calls `mpirun -np N` WITHOUT `--allow-run-as-root` and WITHOUT `--oversubscribe`. The fix is environment variables, not a different OpenMPI version.

**Use system OpenMPI 4.1.6** — do NOT compile 4.1.8 (ABI incompatibility causes `MPI_Type_match_size` error with RIJCOSX).

```bash
ORCA_DIR=/root/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export LD_LIBRARY_PATH=$ORCA_DIR/lib:$ORCA_DIR:$LD_LIBRARY_PATH
export ORCA_TMPDIR=/dev/shm/orca_tmp
export OMPI_MCA_orte_tmpdir_base=/dev/shm/orte_sessions
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_rmaps_base_oversubscribe=1
mkdir -p $ORCA_TMPDIR /dev/shm/orte_sessions
$ORCA_DIR/orca input.inp > output.out 2>&1 &
```

## Input Rules
- `%scf ... end` — MUST close with `end`
- `%pal nprocs N end` — add BEFORE option blocks
- `* xyz charge multiplicity` — odd electrons -> multiplicity 2 (doublet)

## Performance
- 66 atoms PBE0/def2-SVP: ~45s/SCF iteration on 4 cores (vs 160s single-core)
- ~3.5x speedup
