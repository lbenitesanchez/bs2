#!/usr/bin/env Rscript
# Lightweight end-to-end reproducibility check.
# This does not regenerate the full Monte Carlo results from the manuscript.
# It verifies that the code, data, fits, aggregation and plotting interfaces work.

run_cmd <- function(args) {
  cmd <- paste(c("Rscript", args), collapse = " ")
  cat("\n>>> ", cmd, "\n", sep = "")
  status <- system2("Rscript", args = args)
  if (!identical(status, 0L)) stop("Command failed: ", cmd, call. = FALSE)
}

cores <- as.integer(Sys.getenv("CORES", "2"))
backend <- Sys.getenv("BACKEND", "multisession")

cat("BS2 quick reproducibility run\n")
cat("Cores  : ", cores, "\n", sep = "")
cat("Backend: ", backend, "\n", sep = "")

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

run_cmd(c("scripts/00_check_environment.R"))
run_cmd(c("scripts/08_bmd_analysis.R", "--B=19", "--B-sign=999", "--seed=202606", "--out-dir=results/tables_quick"))
run_cmd(c("scripts/01_estimation_parallel.R", "--design=debug", "--R=8", "--chunk=4", paste0("--cores=", cores), paste0("--backend=", backend), "--out-dir=results/raw_estimation_debug", "--prefix=est_debug"))
run_cmd(c("scripts/02_tests_parallel.R", "--design=debug", "--R=8", "--chunk=4", paste0("--cores=", cores), paste0("--backend=", backend), "--B-sign=199", "--out-dir=results/raw_tests_debug", "--prefix=tests_debug"))
run_cmd(c("scripts/03_bootstrap_lr_parallel.R", "--design=debug", "--R=4", "--B=9", "--chunk=2", paste0("--cores=", cores), paste0("--backend=", backend), "--out-dir=results/raw_bootstrap_debug", "--prefix=boot_debug"))
run_cmd(c("scripts/04_aggregate_results.R", "--raw-est=results/raw_estimation_debug", "--raw-tests=results/raw_tests_debug", "--raw-boot=results/raw_bootstrap_debug", "--out-dir=results/tables_debug"))
run_cmd(c("scripts/05_power_parallel.R", "--design=debug", "--R=8", "--chunk=4", paste0("--cores=", cores), paste0("--backend=", backend), "--B-sign=0", "--out=results/raw_power_debug", "--table-dir=results/tables_power_debug"))
run_cmd(c("scripts/06_gof_influence_parallel.R", "--design=debug", "--R=8", "--B=9", "--chunk=2", paste0("--cores=", cores), paste0("--backend=", backend), "--out=results/raw_gof_debug", "--table-dir=results/tables_gof_debug"))
run_cmd(c("scripts/10_algorithmic_benchmark.R", "--R=4", "--chunk=2", paste0("--cores=", cores), paste0("--backend=", backend), "--raw-dir=results/raw_algorithmic_benchmark_debug", "--out-dir=results/tables_algorithmic_benchmark_debug"))

cat("\nQuick run finished. Inspect results/tables_quick, results/tables_debug, results/tables_power_debug and results/tables_gof_debug.\n")
