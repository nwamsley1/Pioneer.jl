#!/bin/bash
#SBATCH --account=compute2-d.goldfarb
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --job-name=lgbm_thread_scaling
set -euo pipefail
BENCH_DIR="${1:?benchmark directory required}"
SOURCE_DIR="${2:?source directory required}"
export JULIA_DEPOT_PATH="$BENCH_DIR/depot:/scratch2/fs1/d.goldfarb/n.t.wamsley/depot"
export JULIA_NUM_THREADS=24 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_NUM_THREADS=64 OMP_DYNAMIC=FALSE
export GKSwstype=100
echo "THREAD BENCHMARK job=$SLURM_JOB_ID node=$SLURMD_NODENAME"
cat "$SOURCE_DIR/PROVENANCE.txt"
lscpu > "$BENCH_DIR/stats/host_cpu.txt"
srun --cpu-bind=cores --export=ALL \
  --container-image=dennisgoldfarb/julia-python:1.12.6 \
  --container-mounts=/storage3/fs1/d.goldfarb/Active:/storage3/fs1/d.goldfarb/Active,/scratch2/fs1/d.goldfarb:/scratch2/fs1/d.goldfarb \
  /usr/local/julia/bin/julia --project="$SOURCE_DIR" --threads=24 --gcthreads=1,1 \
  --compiled-modules=existing "$BENCH_DIR/benchmark.jl"
