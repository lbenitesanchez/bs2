#!/usr/bin/env Rscript
# scripts/04_aggregate_results.R
# -----------------------------------------------------------------------------
# Aggregate raw RDS outputs from the parallel BS2 simulations into CSV tables.
# Run from the project root:
#   Rscript scripts/04_aggregate_results.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source("R/sim_parallel_utils.R")

raw_est_dir <- get_arg("raw-est", "results/raw_estimation")
raw_tests_dir <- get_arg("raw-tests", "results/raw_tests")
raw_boot_dir <- get_arg("raw-boot", "results/raw_bootstrap")
out_dir <- get_arg("out-dir", "results/tables")
ensure_dirs(out_dir)

mcse <- function(p, n) {
  sqrt(pmax(p, 0) * pmax(1 - p, 0) / n)
}

summarise_estimation <- function(raw) {
  if (nrow(raw) == 0L) return(data.table())

  group_cols <- intersect(
    c("experiment", "cell_id", "alpha_id", "n",
      "alpha1_true", "alpha2_true", "beta1_true", "beta2_true", "rho_true"),
    names(raw)
  )

  specs <- list(
    list(method = "ML", prefix = "ml", params = c("alpha1", "alpha2", "beta1", "beta2", "rho"), ok = "ml_ok"),
    list(method = "MM", prefix = "mm", params = c("alpha1", "alpha2", "beta1", "beta2", "rho"), ok = "mm_ok"),
    list(method = "PL", prefix = "pl", params = c("rho"), ok = "pl_ok"),
    list(method = "CF", prefix = "cf", params = c("rho"), ok = "cf_ok")
  )

  out <- list()
  for (sp in specs) {
    for (param in sp$params) {
      est_col <- paste0(sp$prefix, "_", param)
      true_col <- paste0(param, "_true")
      if (!all(c(est_col, true_col) %in% names(raw))) next

      tmp <- raw[
        is.finite(get(est_col)) & is.finite(get(true_col)),
        .(
          R_eff = .N,
          mean_est = mean(get(est_col)),
          true_value = mean(get(true_col)),
          bias = mean(get(est_col) - get(true_col)),
          abs_bias = abs(mean(get(est_col) - get(true_col))),
          rmse = sqrt(mean((get(est_col) - get(true_col))^2)),
          emp_sd = if (.N > 1L) sd(get(est_col)) else NA_real_,
          q025 = as.numeric(quantile(get(est_col), 0.025, na.rm = TRUE, names = FALSE)),
          q500 = as.numeric(quantile(get(est_col), 0.500, na.rm = TRUE, names = FALSE)),
          q975 = as.numeric(quantile(get(est_col), 0.975, na.rm = TRUE, names = FALSE))
        ),
        by = group_cols
      ]

      tmp[, `:=`(
        method = sp$method,
        parameter = param
      )]

      out[[length(out) + 1L]] <- tmp
    }
  }

  if (length(out) == 0L) return(data.table())
  ans <- rbindlist(out, fill = TRUE, use.names = TRUE)
  setcolorder(ans, c(intersect(group_cols, names(ans)), "method", "parameter",
                     setdiff(names(ans), c(group_cols, "method", "parameter"))))
  ans[]
}

summarise_estimation_failures <- function(raw) {
  if (nrow(raw) == 0L) return(data.table())

  group_cols <- intersect(
    c("experiment", "cell_id", "alpha_id", "n",
      "alpha1_true", "alpha2_true", "beta1_true", "beta2_true", "rho_true"),
    names(raw)
  )

  ok_cols <- intersect(c("mm_ok", "ml_ok", "pl_ok", "cf_ok"), names(raw))
  if (length(ok_cols) == 0L) return(data.table())

  out <- lapply(ok_cols, function(ok) {
    raw[, .(
      R_total = .N,
      ok_rate = mean(as.logical(get(ok)), na.rm = TRUE),
      failure_rate = 1 - mean(as.logical(get(ok)), na.rm = TRUE)
    ), by = group_cols][, method := toupper(sub("_ok$", "", ok))]
  })

  rbindlist(out, fill = TRUE, use.names = TRUE)
}

summarise_tests <- function(raw) {
  if (nrow(raw) == 0L) return(data.table())

  group_cols <- intersect(
    c("experiment", "cell_id", "hypothesis", "alt_id", "n", "df", "level",
      "alpha1_true", "alpha2_true", "beta1_true", "beta2_true", "rho_true", "rho0"),
    names(raw)
  )

  reject_cols <- grep("_reject$", names(raw), value = TRUE)
  if (length(reject_cols) == 0L) return(data.table())

  out <- lapply(reject_cols, function(rc) {
    tmp <- raw[!is.na(get(rc)), .(
      R_eff = .N,
      rejection_rate = mean(get(rc)),
      mcse = mcse(mean(get(rc)), .N)
    ), by = group_cols]
    tmp[, statistic := sub("_reject$", "", rc)]
    tmp
  })

  ans <- rbindlist(out, fill = TRUE, use.names = TRUE)
  setcolorder(ans, c(intersect(group_cols, names(ans)), "statistic",
                     setdiff(names(ans), c(group_cols, "statistic"))))
  ans[]
}

summarise_test_fit_rates <- function(raw) {
  if (nrow(raw) == 0L) return(data.table())

  group_cols <- intersect(
    c("experiment", "cell_id", "hypothesis", "alt_id", "n", "df", "level",
      "alpha1_true", "alpha2_true", "beta1_true", "beta2_true", "rho_true", "rho0"),
    names(raw)
  )

  ok_cols <- intersect(c("fit_u_ok", "fit_r_ok", "fit_mm_ok", "tests_ok", "pl_ok"), names(raw))
  if (length(ok_cols) == 0L) return(data.table())

  out <- lapply(ok_cols, function(ok) {
    raw[, .(
      R_total = .N,
      ok_rate = mean(as.logical(get(ok)), na.rm = TRUE),
      failure_rate = 1 - mean(as.logical(get(ok)), na.rm = TRUE)
    ), by = group_cols][, component := sub("_ok$", "", ok)]
  })

  rbindlist(out, fill = TRUE, use.names = TRUE)
}

summarise_bootstrap <- function(raw) {
  if (nrow(raw) == 0L) return(data.table())

  group_cols <- intersect(
    c("experiment", "cell_id", "scenario", "hypothesis", "n", "df", "level",
      "alpha1_true", "alpha2_true", "beta1_true", "beta2_true", "rho_true", "rho0", "B"),
    names(raw)
  )

  raw[, .(
    R_eff = sum(bootstrap_ok == TRUE, na.rm = TRUE),
    R_total = .N,
    bootstrap_ok_rate = mean(bootstrap_ok == TRUE, na.rm = TRUE),
    mean_B_eff = mean(B_eff, na.rm = TRUE),
    mean_inner_failure_rate = mean(bootstrap_failure_rate, na.rm = TRUE),
    asym_rejection_rate = mean(LR_asym_reject, na.rm = TRUE),
    boot_rejection_rate = mean(LR_boot_reject, na.rm = TRUE),
    mcse_asym = mcse(mean(LR_asym_reject, na.rm = TRUE), sum(!is.na(LR_asym_reject))),
    mcse_boot = mcse(mean(LR_boot_reject, na.rm = TRUE), sum(!is.na(LR_boot_reject))),
    mean_LR_obs = mean(LR_obs, na.rm = TRUE),
    mean_boot_crit95 = mean(boot_crit95, na.rm = TRUE),
    median_p_asym = median(p_asym, na.rm = TRUE),
    median_p_boot = median(p_boot, na.rm = TRUE)
  ), by = group_cols]
}

message_header("Aggregating BS2 simulation results")

raw_est <- read_rds_dir(raw_est_dir)
if (nrow(raw_est) > 0L) {
  fwrite(raw_est, file.path(out_dir, "estimation_raw_combined.csv"))
  est_summary <- summarise_estimation(raw_est)
  est_fail <- summarise_estimation_failures(raw_est)
  fwrite(est_summary, file.path(out_dir, "estimation_summary.csv"))
  fwrite(est_fail, file.path(out_dir, "estimation_failure_rates.csv"))
  message("Wrote estimation tables: ", nrow(est_summary), " summary rows.")
} else {
  message("No estimation RDS files found in ", raw_est_dir)
}

raw_tests <- read_rds_dir(raw_tests_dir)
if (nrow(raw_tests) > 0L) {
  fwrite(raw_tests, file.path(out_dir, "tests_raw_combined.csv"))
  test_summary <- summarise_tests(raw_tests)
  test_fit <- summarise_test_fit_rates(raw_tests)
  fwrite(test_summary, file.path(out_dir, "test_rejection_summary.csv"))
  fwrite(test_fit, file.path(out_dir, "test_fit_rates.csv"))
  message("Wrote test tables: ", nrow(test_summary), " rejection summary rows.")
} else {
  message("No test RDS files found in ", raw_tests_dir)
}

raw_boot <- read_rds_dir(raw_boot_dir)
if (nrow(raw_boot) > 0L) {
  fwrite(raw_boot, file.path(out_dir, "bootstrap_raw_combined.csv"))
  boot_summary <- summarise_bootstrap(raw_boot)
  fwrite(boot_summary, file.path(out_dir, "bootstrap_summary.csv"))
  message("Wrote bootstrap tables: ", nrow(boot_summary), " summary rows.")
} else {
  message("No bootstrap RDS files found in ", raw_boot_dir)
}

message("Aggregation finished. Tables are in: ", out_dir)
