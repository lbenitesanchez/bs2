#!/usr/bin/env bash
set -euo pipefail

# Full Ubuntu reproduction workflow for the BS2 manuscript computations.
# Run from the repository root, preferably inside tmux/screen/nohup.
# The full workflow can take many hours depending on cores and CPU speed.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

CORES=${CORES:-$(( $(nproc) - 1 ))}
if [ "${CORES}" -lt 1 ]; then CORES=1; fi
BACKEND=${BACKEND:-multicore}

mkdir -p logs results figures

echo "Full BS2 workflow"
echo "CORES=${CORES}"
echo "BACKEND=${BACKEND}"
echo "Started at $(date -Is)"

Rscript scripts/00_check_environment.R | tee logs/00_check_environment.log

# Empirical application tables.
Rscript scripts/08_bmd_analysis.R \
  --B=499 \
  --B-sign=99999 \
  --seed=202606 \
  --out-dir=results/tables \
  > logs/08_bmd_analysis.log 2>&1

# Monte Carlo experiment 1: point estimation.
Rscript scripts/01_estimation_parallel.R \
  --design=estimation \
  --R=1500 \
  --chunk=25 \
  --cores=${CORES} \
  --backend=${BACKEND} \
  --seed=202501 \
  --out-dir=results/raw_estimation \
  --prefix=estimation \
  > logs/01_estimation_parallel.log 2>&1

# Monte Carlo experiment 2: size.
Rscript scripts/02_tests_parallel.R \
  --design=size \
  --R=3000 \
  --chunk=25 \
  --cores=${CORES} \
  --backend=${BACKEND} \
  --B-sign=9999 \
  --seed=202502 \
  --out-dir=results/raw_tests \
  --prefix=tests_size \
  > logs/02_tests_size_parallel.log 2>&1

# Bootstrap LR stress experiment.
Rscript scripts/03_bootstrap_lr_parallel.R \
  --design=stress \
  --R=500 \
  --B=499 \
  --chunk=5 \
  --cores=${CORES} \
  --backend=${BACKEND} \
  --seed=202503 \
  --out-dir=results/raw_bootstrap \
  --prefix=boot_stress \
  > logs/03_bootstrap_lr_parallel.log 2>&1

# Power experiment.
Rscript scripts/05_power_parallel.R \
  --design=power \
  --R=1500 \
  --chunk=25 \
  --cores=${CORES} \
  --backend=${BACKEND} \
  --B-sign=9999 \
  --seed=202504 \
  --out=results/raw_power_score_fixed \
  --table-dir=results/tables_power_score_fixed \
  > logs/05_power_parallel.log 2>&1

# GOF/influence simulation.
Rscript scripts/06_gof_influence_parallel.R \
  --design=gof \
  --R=500 \
  --B=199 \
  --chunk=5 \
  --cores=${CORES} \
  --backend=${BACKEND} \
  --seed=202505 \
  --out=results/raw_gof_influence \
  --table-dir=results/tables_gof_influence \
  > logs/06_gof_influence_parallel.log 2>&1


# Algorithmic benchmark: profiled 2D optimizer versus direct 5D optimizer.
Rscript scripts/10_algorithmic_benchmark.R \
  --R=200 \
  --chunk=10 \
  --cores=${CORES} \
  --backend=${BACKEND} \
  --seed=202507 \
  --raw-dir=results/raw_algorithmic_benchmark \
  --out-dir=results/tables_algorithmic_benchmark \
  > logs/10_algorithmic_benchmark.log 2>&1

# Aggregate raw RDS files produced by scripts 01--03.
Rscript scripts/04_aggregate_results.R \
  --raw-est=results/raw_estimation \
  --raw-tests=results/raw_tests \
  --raw-boot=results/raw_bootstrap \
  --out-dir=results/tables \
  > logs/04_aggregate_results.log 2>&1

# Publication figures. Use --formats=eps,pdf for local preview.
Rscript scripts/07_make_figures_publication.R \
  --fig-dir=figures \
  --formats=eps,pdf \
  --B-env=500 \
  --seed=202506 \
  --bmd-file=data/bmd.csv \
  > logs/07_make_figures_publication.log 2>&1

Rscript scripts/00_check_environment.R > logs/final_session_info.log 2>&1

echo "Finished at $(date -Is)"
