# BS2 reproducibility bundle

This repository contains the R code and data used for the paper:

**Computationally stable inference and diagnostics for the bivariate Birnbaum--Saunders distribution**

The code implements profiled maximum likelihood, modified-moment estimation, conditional sign pseudo-likelihood, LR/score/Wald/MM-Wald tests, parametric-bootstrap LR calibration, BMD diagnostics, and the Monte Carlo workflows reported in the manuscript and supplement.

## Repository structure

```text
bs2/
  R/
    bs2_core.R                 # BS2 density/log-likelihood, simulation, MM/ML/PL fits, tests, bootstrap helpers
    sim_designs.R              # Monte Carlo design grids
    sim_parallel_utils.R       # CLI parsing, reproducible seeds, future backends, RDS checkpointing
  data/
    bmd.csv                    # BMD data used in the empirical application
  scripts/
    00_install.R               # install R package dependencies
    00_check_environment.R     # check packages/parallel support and write sessionInfo
    00_freeze_renv.R           # optional renv.lock creation on the final analysis machine
    00_run_all_quick.R         # fast smoke/end-to-end reproducibility check
    00_run_all_full_ubuntu.sh  # full Ubuntu workflow for manuscript-scale computations
    01_estimation_parallel.R   # Monte Carlo point-estimation experiment
    02_tests_parallel.R        # Monte Carlo size experiment
    03_bootstrap_lr_parallel.R # LR bootstrap stress experiment
    04_aggregate_results.R     # aggregate raw RDS outputs into CSV tables
    05_power_parallel.R        # Monte Carlo power experiment
    06_gof_influence_parallel.R# GOF and influence simulation
    07_make_figures_publication.R
    08_bmd_analysis.R          # empirical BMD tables
    09_bmd_bs2_lognormal_contours.R
    10_algorithmic_benchmark.R # profile-vs-full optimizer benchmark
  manuscript/
    table_figure_map.csv       # mapping from manuscript outputs to scripts
  results/                     # generated outputs; raw files are ignored by git
  figures/                     # generated figures; ignored by git except README/.gitkeep
  logs/                        # run logs; ignored by git except README/.gitkeep
  tests/
    smoke_test.R               # short numerical sanity check
```

## System requirements

The code was organized for Ubuntu/Linux parallel execution. It also runs sequentially or with multisession workers on other systems, but `--backend=multicore` is intended for Ubuntu/Linux.

Recommended Ubuntu packages:

```bash
sudo apt update
sudo apt install -y r-base r-base-dev git make \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libcairo2-dev libxt-dev libfontconfig1-dev
```

Required R packages:

```r
c("data.table", "future", "future.apply", "parallelly", "RhpcBLASctl", "ggplot2")
```

Install them with:

```bash
Rscript scripts/00_install.R
```

Alternatively, create a conda environment from `environment.yml` if your workflow uses conda/mamba.

## Quick reproducibility check

From the repository root:

```bash
Rscript scripts/00_check_environment.R
CORES=2 BACKEND=multisession Rscript scripts/00_run_all_quick.R
Rscript tests/smoke_test.R
```

The quick workflow uses tiny debug grids. It does **not** reproduce the manuscript Monte Carlo tables, but it verifies that the repository, data, fits, checkpointed parallel tasks, aggregation and BMD analysis pipeline work on a fresh machine.

## Full Ubuntu reproduction

Run inside `tmux`, `screen` or `nohup` because the full Monte Carlo/bootstrapping workflow can be long.

```bash
chmod +x scripts/00_run_all_full_ubuntu.sh
CORES=$(($(nproc)-1)) BACKEND=multicore ./scripts/00_run_all_full_ubuntu.sh
```

The full script runs:

1. BMD empirical analysis and bootstrap LR p-values.
2. Point-estimation Monte Carlo experiment.
3. Size Monte Carlo experiment.
4. Bootstrap LR stress experiment.
5. Power Monte Carlo experiment.
6. GOF/influence Monte Carlo experiment.
7. Profile-vs-full optimizer algorithmic benchmark.
8. Aggregation of raw RDS files to CSV.
9. Publication figures in EPS/PDF.
10. Final `sessionInfo()` capture.

The raw Monte Carlo files are checkpointed under `results/raw_*`. If a run is interrupted, re-running the same command resumes only pending chunks because completed chunk files are detected automatically.

## Parallel and random-number controls

The scripts use `future.apply` with either `multicore`, `multisession` or `sequential` backends. The helper `R/sim_parallel_utils.R` sets `RNGkind("L'Ecuyer-CMRG")`, creates deterministic cell/chunk seeds, and pins BLAS/OpenMP numerical threads to one thread per worker to avoid nested parallelism.

Typical Ubuntu command:

```bash
Rscript scripts/01_estimation_parallel.R \
  --design=estimation \
  --R=1500 \
  --chunk=25 \
  --cores=14 \
  --backend=multicore \
  --seed=202501 \
  --out-dir=results/raw_estimation \
  --prefix=estimation
```

For laptops, GitHub Actions, Windows or macOS, prefer:

```bash
--backend=multisession
```

For strict single-process debugging:

```bash
--cores=1 --backend=sequential
```

## Reproducing manuscript outputs

The file `manuscript/table_figure_map.csv` maps each main/supplementary table or figure to the script(s) that generate the corresponding raw outputs, summaries or figures.

Common commands:

```bash
# Aggregate the main raw outputs created by scripts 01--03
Rscript scripts/04_aggregate_results.R \
  --raw-est=results/raw_estimation \
  --raw-tests=results/raw_tests \
  --raw-boot=results/raw_bootstrap \
  --out-dir=results/tables

# Empirical BMD tables
Rscript scripts/08_bmd_analysis.R --B=499 --B-sign=99999 --out-dir=results/tables

# Publication figures
Rscript scripts/07_make_figures_publication.R --fig-dir=figures --formats=eps,pdf

# Supplementary BS2/lognormal BMD contours
Rscript scripts/09_bmd_bs2_lognormal_contours.R --fig-dir=figures --formats=pdf,png

# Algorithmic profile-vs-full optimizer benchmark
Rscript scripts/10_algorithmic_benchmark.R --R=200 --cores=14 --backend=multicore
```

## Version locking

This repository includes a `DESCRIPTION` file listing all R dependencies. To freeze the exact package versions on the final analysis machine, run:

```bash
Rscript scripts/00_freeze_renv.R
```

This creates or updates `renv.lock`. Commit the resulting lockfile before journal submission if exact package-version replay is required.

## Git and large files

The `.gitignore` keeps raw Monte Carlo chunks, logs and generated figures out of git by default. For peer review, commit code, data, documentation, small CSV summary tables if desired, and `renv.lock` after freezing. Large raw simulation outputs should be archived as a release asset, Zenodo deposit, OSF component or similar reproducibility archive.

## License

Code is released under the MIT License. The BMD data are documented separately in `data/README.md`; the code license does not automatically license the data.
