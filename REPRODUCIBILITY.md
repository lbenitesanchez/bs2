# Reproducibility guide

This project is designed as a checkpointed, parallel R workflow for Ubuntu/Linux.

## Minimal validation

Use this before pushing to GitHub or creating a release:

```bash
Rscript scripts/00_install.R
Rscript scripts/00_check_environment.R
CORES=2 BACKEND=multisession Rscript scripts/00_run_all_quick.R
Rscript tests/smoke_test.R
```

The quick run writes small debug outputs only. It checks that the computational core can be sourced, the BMD data can be read, the unrestricted/restricted fits run, the bootstrap wrapper runs on a tiny B, and the Monte Carlo scripts can checkpoint/aggregate.

## Full manuscript-scale validation

```bash
CORES=$(($(nproc)-1)) BACKEND=multicore ./scripts/00_run_all_full_ubuntu.sh
```

The full workflow uses the replication counts reported in the manuscript/supplement:

- point estimation: `R = 1500` per design cell;
- empirical size: `R = 3000` per design cell;
- empirical power: `R = 1500` per design cell;
- LR bootstrap stress: outer `R = 500`, inner `B = 499`;
- GOF/influence simulation: outer `R = 500`, inner `B = 199`;
- algorithmic benchmark: `R = 200` per benchmark cell;
- BMD LR bootstrap: `B = 499`;
- BMD sign resampling: `B = 99999`.

## Seeds and parallelism

The simulation workflow defines tasks by design cell and replication block. Each task has a deterministic seed generated from the base seed, design cell and block index. Completed task files are skipped when rerunning the same command, so interrupted runs can be resumed.

`R/sim_parallel_utils.R` sets:

```r
RNGkind("L'Ecuyer-CMRG")
```

and pins numerical threads to one thread per worker through environment variables and `RhpcBLASctl` when available. This avoids oversubscription when many R workers are running.

## Recommended release contents

For journal submission or archival release, include:

1. all source code in `R/` and `scripts/`;
2. `data/bmd.csv` and `data/README.md`;
3. `README.md`, `REPRODUCIBILITY.md`, `DESCRIPTION`, `CITATION.cff`, `LICENSE`;
4. `renv.lock`, generated on the final analysis machine with `scripts/00_freeze_renv.R`;
5. `results/sessionInfo.txt` generated after the final run;
6. small CSV summaries needed by the manuscript, if the repository will not include raw RDS chunks;
7. a GitHub release or external archive for large raw Monte Carlo outputs, if exact raw replay is desired without recomputation.

## Why raw outputs are ignored by default

Full raw Monte Carlo and bootstrap outputs can be large and are expensive to regenerate. The `.gitignore` excludes raw RDS chunks and generated figures/tables by default. Commit final summary CSVs selectively, or archive large results in a release asset or external repository.
