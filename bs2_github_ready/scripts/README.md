# Scripts

Run scripts from the repository root. Most scripts accept command-line arguments of the form `--name=value`.

Use `00_run_all_quick.R` for a short validation run and `00_run_all_full_ubuntu.sh` for the manuscript-scale Ubuntu workflow. The full workflow is parallelized and checkpointed; reruns skip already completed RDS chunks.

Main manuscript-scale scripts:

- `01_estimation_parallel.R`: point-estimation Monte Carlo.
- `02_tests_parallel.R`: empirical size Monte Carlo.
- `03_bootstrap_lr_parallel.R`: LR bootstrap calibration stress experiment.
- `04_aggregate_results.R`: aggregation of raw RDS chunks from scripts 01--03.
- `05_power_parallel.R`: empirical power Monte Carlo.
- `06_gof_influence_parallel.R`: bootstrap GOF and influence simulation.
- `07_make_figures_publication.R`: publication figures from data and aggregated tables.
- `08_bmd_analysis.R`: empirical BMD estimates, tests, bootstrap and model comparison tables.
- `09_bmd_bs2_lognormal_contours.R`: supplementary BMD contour comparison.
- `10_algorithmic_benchmark.R`: profile-likelihood versus full 5D optimization benchmark.
