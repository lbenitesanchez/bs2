#!/usr/bin/env Rscript
# scripts/06_gof_influence_parallel.R
# -----------------------------------------------------------------------------
# BS2 goodness-of-fit and influence/sensitivity simulation.
#
# This is the E5 diagnostic simulation complementing the estimation, size,
# power and bootstrap-LR experiments.  It studies the bootstrap Cramer--von
# Mises goodness-of-fit statistic based on fitted Mahalanobis distances and
# a case-deletion sensitivity summary under:
#   1) correctly specified BS2 data,
#   2) bivariate lognormal data, and
#   3) contaminated BS2 data with a small low-scale component.
#
# Typical smoke test:
#   Rscript scripts/06_gof_influence_parallel.R \
#     --design=debug --R=2 --B=9 --chunk=1 --cores=1 --backend=sequential
#
# Typical formal run:
#   CORES=$(( $(nproc) - 2 )); if [ "$CORES" -lt 1 ]; then CORES=1; fi
#   nohup Rscript scripts/06_gof_influence_parallel.R \
#     --design=gof --R=500 --B=199 --chunk=5 --cores=${CORES} \
#     --backend=multicore --seed=202505 \
#     --out=results/raw_gof_influence \
#     --table-dir=results/tables_gof_influence \
#     > logs/gof_influence.log 2>&1 &
#
# If compute budget permits, increase --B=499 for the final run.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

source("R/bs2_core.R")
source("R/sim_parallel_utils.R")
source("R/sim_designs.R")

# -----------------------------------------------------------------------------
# Compatibility fallbacks for older project folders
# -----------------------------------------------------------------------------

if (!exists("log_message", mode = "function")) {
  log_message <- function(...) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
  }
}

if (!exists("write_rds_atomic", mode = "function")) {
  write_rds_atomic <- function(object, file) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(file, ".tmp_", Sys.getpid(), "_", sample.int(1e9, 1L))
    saveRDS(object, tmp)
    ok <- file.rename(tmp, file)
    if (!ok) {
      file.copy(tmp, file, overwrite = TRUE)
      unlink(tmp)
    }
    invisible(file)
  }
}

if (!exists("read_rds_dir", mode = "function")) {
  read_rds_dir <- function(path, pattern = "\\.rds$") {
    files <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0L) return(data.table())
    rbindlist(lapply(files, readRDS), fill = TRUE)
  }
}

if (!exists("get_bool_arg", mode = "function")) {
  get_bool_arg <- function(name, default = FALSE) {
    val <- tolower(get_arg(name, if (isTRUE(default)) "true" else "false"))
    if (val %in% c("true", "t", "1", "yes", "y")) return(TRUE)
    if (val %in% c("false", "f", "0", "no", "n")) return(FALSE)
    stop("Argument --", name, " must be TRUE/FALSE, yes/no, or 1/0.")
  }
}

parse_csv_character <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_csv_integer <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  out <- suppressWarnings(as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1L]])))
  if (any(is.na(out))) stop("Integer CSV argument contains NA: ", x)
  out
}

# -----------------------------------------------------------------------------
# Command-line arguments
# -----------------------------------------------------------------------------

cores <- get_int_arg("cores", max(1L, parallelly::availableCores() - 1L))
R <- get_int_arg("R", 500L)
B <- get_int_arg("B", 199L)
chunk_size <- get_int_arg("chunk", 5L)
seed <- get_int_arg("seed", 202505L)
backend <- get_arg("backend", "multicore")
design_name <- tolower(get_arg("design", "gof"))
debug <- get_bool_arg("debug", FALSE) || identical(design_name, "debug")
aggregate_after <- get_bool_arg("aggregate", TRUE)
do_influence <- get_bool_arg("influence", TRUE)
do_gof_boot <- get_bool_arg("gof-bootstrap", TRUE)

n_grid_arg <- get_arg("n-grid", NULL)
dgp_arg <- get_arg("dgp", NULL)
n_grid <- parse_csv_integer(n_grid_arg, default = c(24L, 50L, 100L))
dgp_grid <- parse_csv_character(dgp_arg, default = c("bs2", "lognormal", "contaminated_bs2"))

out_arg <- get_arg("out", get_arg("out-dir", NULL))
out_dir <- if (!is.null(out_arg)) {
  out_arg
} else if (debug) {
  "results/raw_gof_debug"
} else {
  "results/raw_gof_influence"
}

table_arg <- get_arg("table-dir", NULL)
table_dir <- if (!is.null(table_arg)) {
  table_arg
} else if (debug) {
  "results/tables_gof_debug"
} else {
  "results/tables_gof_influence"
}

prefix <- get_arg("prefix", if (debug) "gof_debug" else "gof")

setup_parallel(cores = cores, backend = backend)

# -----------------------------------------------------------------------------
# Design construction
# -----------------------------------------------------------------------------

if (debug) {
  design <- make_gof_design(n_grid = c(24L), dgp_grid = c("bs2", "contaminated_bs2"))
  R <- min(R, 8L)
  B <- min(B, 9L)
  chunk_size <- min(chunk_size, 2L)
} else {
  if (!design_name %in% c("gof", "e5", "diagnostic", "gof_influence")) {
    stop("Unknown --design value: ", design_name, ". Use gof or debug.")
  }
  design <- make_gof_design(n_grid = n_grid, dgp_grid = dgp_grid)
}

design <- add_cell_id(design)

tasks <- make_rep_tasks(
  design = design,
  R = R,
  chunk_size = chunk_size,
  out_dir = out_dir,
  prefix = prefix
)

tasks_todo <- tasks[!file.exists(out_file)]

log_message("GOF/influence design cells: ", nrow(design))
log_message("Total tasks: ", nrow(tasks))
log_message("Pending tasks: ", nrow(tasks_todo))
log_message("R per cell: ", R)
log_message("Bootstrap B per replication: ", B)
log_message("Chunk size: ", chunk_size)
log_message("Workers: ", cores)
log_message("Backend: ", backend)
log_message("Do GOF bootstrap: ", do_gof_boot)
log_message("Do influence: ", do_influence)
log_message("Output directory: ", out_dir)
log_message("Table directory: ", table_dir)

# -----------------------------------------------------------------------------
# Statistical helpers
# -----------------------------------------------------------------------------

bs2_mahalanobis <- function(t1, t2, alpha, beta, rho) {
  z1 <- a_bs(t1, alpha[1L], beta[1L])
  z2 <- a_bs(t2, alpha[2L], beta[2L])
  (z1^2 + z2^2 - 2 * rho * z1 * z2) / (1 - rho^2)
}

cvm_chisq2_stat <- function(d) {
  d <- sort(as.numeric(d))
  n <- length(d)
  u <- pchisq(d, df = 2)
  1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n))^2)
}

fit_bs2_full <- function(t1, t2) {
  fit_mm <- safe_eval(fit_bs2_mm(t1, t2))
  if (!fit_ok(fit_mm)) {
    return(structure(list(error = "MM fit failed"), class = "try-error"))
  }
  fit_ml <- safe_eval(fit_bs2_profile(t1, t2, beta_start = fit_mm$beta))
  if (!fit_ok(fit_ml)) {
    return(structure(list(error = "ML fit failed"), class = "try-error"))
  }
  fit_ml$fit_mm <- fit_mm
  fit_ml
}

compute_gof_stat <- function(t1, t2, fit_ml) {
  d <- bs2_mahalanobis(t1, t2, fit_ml$alpha, fit_ml$beta, fit_ml$rho)
  list(
    d = d,
    W2 = cvm_chisq2_stat(d),
    max_d = max(d),
    median_d = median(d)
  )
}

# Bivariate lognormal alternative calibrated to the BMD data on the log scale.
.bmd_for_lognormal <- make_bmd_data()
.log_bmd <- as.matrix(log(.bmd_for_lognormal[, .(t1, t2)]))
.lognorm_mu <- colMeans(.log_bmd)
.lognorm_Sigma <- cov(.log_bmd)
# Numerical guard against a nearly singular covariance estimate.
.lognorm_Sigma <- .lognorm_Sigma + diag(1e-10, 2L)

rlognormal_biv_bmd <- function(n) {
  Rchol <- chol(.lognorm_Sigma)
  z <- matrix(rnorm(2L * n), nrow = n, ncol = 2L) %*% Rchol
  z <- sweep(z, 2L, .lognorm_mu, FUN = "+")
  exp(z)
}

simulate_gof_dgp <- function(n, dgp, alpha, beta, rho, eps_contam, beta_contam) {
  dgp <- as.character(dgp)
  if (dgp == "bs2") {
    return(rbs2(n = n, alpha = alpha, beta = beta, rho = rho))
  }

  if (dgp == "lognormal") {
    return(rlognormal_biv_bmd(n))
  }

  if (dgp == "contaminated_bs2") {
    base <- rbs2(n = n, alpha = alpha, beta = beta, rho = rho)
    contam <- rbs2(n = n, alpha = alpha, beta = beta_contam, rho = rho)
    is_contam <- runif(n) < eps_contam
    base[is_contam, ] <- contam[is_contam, ]
    return(base)
  }

  stop("Unknown DGP: ", dgp)
}

gof_bootstrap_refit <- function(t1, t2, fit_ml, W2_obs, B) {
  if (B <= 0L) {
    return(list(p_boot = NA_real_, crit95 = NA_real_, B_eff = 0L,
                boot_fail_rate = NA_real_))
  }

  n <- length(t1)
  W2_star <- rep(NA_real_, B)
  ok <- rep(FALSE, B)

  for (b in seq_len(B)) {
    dat_b <- rbs2(n = n, alpha = fit_ml$alpha, beta = fit_ml$beta, rho = fit_ml$rho)
    fit_b <- fit_bs2_full(dat_b[, 1L], dat_b[, 2L])
    if (fit_ok(fit_b)) {
      stat_b <- safe_eval(compute_gof_stat(dat_b[, 1L], dat_b[, 2L], fit_b))
      if (fit_ok(stat_b) && is.finite(stat_b$W2)) {
        W2_star[b] <- stat_b$W2
        ok[b] <- TRUE
      }
    }
  }

  B_eff <- sum(ok)
  if (B_eff == 0L) {
    return(list(p_boot = NA_real_, crit95 = NA_real_, B_eff = 0L,
                boot_fail_rate = 1.0))
  }

  W2_ok <- W2_star[ok]
  list(
    p_boot = (1 + sum(W2_ok >= W2_obs)) / (B_eff + 1),
    crit95 = as.numeric(quantile(W2_ok, 0.95, type = 1, names = FALSE)),
    B_eff = B_eff,
    boot_fail_rate = 1 - B_eff / B
  )
}

case_deletion_summary <- function(t1, t2, fit_full) {
  n <- length(t1)
  theta_full <- c(fit_full$alpha, fit_full$beta, fit_full$rho)
  denom <- pmax(abs(theta_full), 1e-8)

  rel_norm <- rep(NA_real_, n)
  delta_alpha1 <- delta_alpha2 <- delta_beta1 <- delta_beta2 <- delta_rho <- rep(NA_real_, n)
  ok <- rep(FALSE, n)

  for (i in seq_len(n)) {
    idx <- setdiff(seq_len(n), i)
    fit_i <- fit_bs2_full(t1[idx], t2[idx])
    if (fit_ok(fit_i)) {
      theta_i <- c(fit_i$alpha, fit_i$beta, fit_i$rho)
      delta <- theta_i - theta_full
      rel_norm[i] <- sqrt(sum((delta / denom)^2))
      delta_alpha1[i] <- delta[1L]
      delta_alpha2[i] <- delta[2L]
      delta_beta1[i] <- delta[3L]
      delta_beta2[i] <- delta[4L]
      delta_rho[i] <- delta[5L]
      ok[i] <- TRUE
    }
  }

  if (!any(ok)) {
    return(list(
      influence_ok = FALSE,
      max_rel_norm = NA_real_, median_rel_norm = NA_real_, mean_rel_norm = NA_real_,
      top_case = NA_integer_, second_rel_norm = NA_real_, deletion_fail_rate = 1.0,
      top_delta_alpha1 = NA_real_, top_delta_alpha2 = NA_real_,
      top_delta_beta1 = NA_real_, top_delta_beta2 = NA_real_, top_delta_rho = NA_real_
    ))
  }

  ord <- order(rel_norm, decreasing = TRUE, na.last = NA)
  top <- ord[1L]
  second <- if (length(ord) >= 2L) rel_norm[ord[2L]] else NA_real_

  list(
    influence_ok = TRUE,
    max_rel_norm = rel_norm[top],
    median_rel_norm = median(rel_norm[ok], na.rm = TRUE),
    mean_rel_norm = mean(rel_norm[ok], na.rm = TRUE),
    top_case = top,
    second_rel_norm = second,
    deletion_fail_rate = mean(!ok),
    top_delta_alpha1 = delta_alpha1[top],
    top_delta_alpha2 = delta_alpha2[top],
    top_delta_beta1 = delta_beta1[top],
    top_delta_beta2 = delta_beta2[top],
    top_delta_rho = delta_rho[top]
  )
}

one_gof_replication <- function(cell, rep_id) {
  pin_numeric_threads()

  n <- as.integer(cell$n)
  dgp <- as.character(cell$dgp)
  alpha <- c(as.numeric(cell$alpha1), as.numeric(cell$alpha2))
  beta <- c(as.numeric(cell$beta1), as.numeric(cell$beta2))
  rho <- as.numeric(cell$rho_true)
  eps_contam <- as.numeric(cell$eps_contam)
  beta_contam <- c(as.numeric(cell$beta1_contam), as.numeric(cell$beta2_contam))

  dat <- simulate_gof_dgp(
    n = n, dgp = dgp, alpha = alpha, beta = beta, rho = rho,
    eps_contam = eps_contam, beta_contam = beta_contam
  )
  t1 <- dat[, 1L]
  t2 <- dat[, 2L]

  fit_ml <- fit_bs2_full(t1, t2)
  ml_ok <- fit_ok(fit_ml)

  W2 <- max_d <- median_d <- NA_real_
  p_boot <- crit95 <- NA_real_
  B_eff <- 0L
  boot_fail_rate <- NA_real_
  gof_reject <- NA_integer_

  infl <- list(
    influence_ok = NA, max_rel_norm = NA_real_, median_rel_norm = NA_real_,
    mean_rel_norm = NA_real_, top_case = NA_integer_, second_rel_norm = NA_real_,
    deletion_fail_rate = NA_real_, top_delta_alpha1 = NA_real_, top_delta_alpha2 = NA_real_,
    top_delta_beta1 = NA_real_, top_delta_beta2 = NA_real_, top_delta_rho = NA_real_
  )

  if (ml_ok) {
    stat <- safe_eval(compute_gof_stat(t1, t2, fit_ml))
    if (fit_ok(stat)) {
      W2 <- stat$W2
      max_d <- stat$max_d
      median_d <- stat$median_d
    }

    if (isTRUE(do_gof_boot) && is.finite(W2)) {
      boot <- safe_eval(gof_bootstrap_refit(t1, t2, fit_ml, W2_obs = W2, B = B))
      if (fit_ok(boot)) {
        p_boot <- boot$p_boot
        crit95 <- boot$crit95
        B_eff <- boot$B_eff
        boot_fail_rate <- boot$boot_fail_rate
        gof_reject <- if (is.finite(p_boot)) as.integer(p_boot < 0.05) else NA_integer_
      }
    }

    if (isTRUE(do_influence)) {
      infl_tmp <- safe_eval(case_deletion_summary(t1, t2, fit_ml))
      if (fit_ok(infl_tmp)) infl <- infl_tmp
    }
  }

  data.table(
    rep = rep_id,
    experiment = ifelse("experiment" %in% names(cell), cell$experiment, "E5_gof_influence"),
    cell_id = ifelse("cell_id" %in% names(cell), cell$cell_id, NA_integer_),
    dgp = dgp,
    model_correct = dgp == "bs2",
    n = n,
    alpha1_true = alpha[1L],
    alpha2_true = alpha[2L],
    beta1_true = beta[1L],
    beta2_true = beta[2L],
    rho_true = rho,
    eps_contam = eps_contam,
    beta1_contam = beta_contam[1L],
    beta2_contam = beta_contam[2L],

    ml_ok = ml_ok,
    alpha1_hat = if (ml_ok) fit_ml$alpha[1L] else NA_real_,
    alpha2_hat = if (ml_ok) fit_ml$alpha[2L] else NA_real_,
    beta1_hat = if (ml_ok) fit_ml$beta[1L] else NA_real_,
    beta2_hat = if (ml_ok) fit_ml$beta[2L] else NA_real_,
    rho_hat = if (ml_ok) fit_ml$rho else NA_real_,
    logLik = if (ml_ok) fit_ml$logLik else NA_real_,

    W2 = W2,
    max_d = max_d,
    median_d = median_d,
    B = B,
    B_eff = B_eff,
    p_boot = p_boot,
    crit95 = crit95,
    boot_fail_rate = boot_fail_rate,
    gof_reject = gof_reject,

    influence_ok = as.logical(infl$influence_ok),
    max_rel_norm = as.numeric(infl$max_rel_norm),
    median_rel_norm = as.numeric(infl$median_rel_norm),
    mean_rel_norm = as.numeric(infl$mean_rel_norm),
    top_case = as.integer(infl$top_case),
    second_rel_norm = as.numeric(infl$second_rel_norm),
    deletion_fail_rate = as.numeric(infl$deletion_fail_rate),
    top_delta_alpha1 = as.numeric(infl$top_delta_alpha1),
    top_delta_alpha2 = as.numeric(infl$top_delta_alpha2),
    top_delta_beta1 = as.numeric(infl$top_delta_beta1),
    top_delta_beta2 = as.numeric(infl$top_delta_beta2),
    top_delta_rho = as.numeric(infl$top_delta_rho)
  )
}

run_gof_task <- function(task_row) {
  if (file.exists(task_row$out_file)) return(task_row$out_file)

  cell <- design[task_row$cell_id]
  set.seed(task_seed(seed, task_row$cell_id, task_row$chunk_id))
  reps <- seq.int(task_row$rep_start, task_row$rep_end)

  ans <- lapply(reps, function(r) one_gof_replication(cell = cell, rep_id = r))
  ans <- rbindlist(ans, fill = TRUE)
  write_rds_atomic(ans, task_row$out_file)
  task_row$out_file
}

# -----------------------------------------------------------------------------
# Run checkpointed parallel tasks
# -----------------------------------------------------------------------------

if (nrow(tasks_todo) > 0L) {
  invisible(future_lapply(
    seq_len(nrow(tasks_todo)),
    function(ii) run_gof_task(tasks_todo[ii]),
    future.seed = TRUE
  ))
}

log_message("GOF/influence simulation finished.")

# -----------------------------------------------------------------------------
# Aggregation
# -----------------------------------------------------------------------------

aggregate_gof_results <- function(raw_dir, table_dir) {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  raw <- read_rds_dir(raw_dir)

  if (nrow(raw) == 0L) {
    warning("No RDS files found in ", raw_dir, ". No tables were written.")
    return(invisible(NULL))
  }

  group_cols <- intersect(
    c("experiment", "dgp", "model_correct", "n", "alpha1_true", "alpha2_true",
      "beta1_true", "beta2_true", "rho_true", "eps_contam", "beta1_contam", "beta2_contam"),
    names(raw)
  )

  gof_summary <- raw[, .(
    R_total = .N,
    ml_success_rate = mean(ml_ok, na.rm = TRUE),
    mean_W2 = mean(W2, na.rm = TRUE),
    median_W2 = median(W2, na.rm = TRUE),
    mean_max_d = mean(max_d, na.rm = TRUE),
    median_max_d = median(max_d, na.rm = TRUE),
    mean_B_eff = mean(B_eff, na.rm = TRUE),
    mean_boot_fail_rate = mean(boot_fail_rate, na.rm = TRUE),
    rejection_rate = mean(gof_reject, na.rm = TRUE),
    rejection_mcse = sqrt(mean(gof_reject, na.rm = TRUE) * (1 - mean(gof_reject, na.rm = TRUE)) / sum(!is.na(gof_reject))),
    median_p_boot = median(p_boot, na.rm = TRUE),
    mean_crit95 = mean(crit95, na.rm = TRUE),
    R_gof_eff = sum(!is.na(gof_reject))
  ), by = group_cols]

  setorder(gof_summary, dgp, n)
  fwrite(gof_summary, file.path(table_dir, "gof_rejection_summary.csv"))

  influence_summary <- raw[, .(
    R_total = .N,
    influence_success_rate = mean(influence_ok, na.rm = TRUE),
    mean_max_rel_norm = mean(max_rel_norm, na.rm = TRUE),
    median_max_rel_norm = median(max_rel_norm, na.rm = TRUE),
    q90_max_rel_norm = as.numeric(quantile(max_rel_norm, 0.90, na.rm = TRUE, names = FALSE)),
    q95_max_rel_norm = as.numeric(quantile(max_rel_norm, 0.95, na.rm = TRUE, names = FALSE)),
    mean_second_rel_norm = mean(second_rel_norm, na.rm = TRUE),
    median_second_rel_norm = median(second_rel_norm, na.rm = TRUE),
    mean_deletion_fail_rate = mean(deletion_fail_rate, na.rm = TRUE),
    mean_abs_top_delta_rho = mean(abs(top_delta_rho), na.rm = TRUE),
    median_abs_top_delta_rho = median(abs(top_delta_rho), na.rm = TRUE)
  ), by = group_cols]

  setorder(influence_summary, dgp, n)
  fwrite(influence_summary, file.path(table_dir, "influence_summary.csv"))

  fit_summary <- raw[, .(
    R_total = .N,
    ml_success_rate = mean(ml_ok, na.rm = TRUE),
    ml_failure_rate = mean(!ml_ok, na.rm = TRUE),
    gof_effective_rate = mean(!is.na(gof_reject)),
    influence_effective_rate = mean(!is.na(influence_ok))
  ), by = group_cols]

  setorder(fit_summary, dgp, n)
  fwrite(fit_summary, file.path(table_dir, "gof_fit_rates.csv"))

  # Wide convenience table for the manuscript.
  manuscript <- merge(
    gof_summary[, .(dgp, model_correct, n, R_gof_eff, rejection_rate, rejection_mcse, median_p_boot, mean_W2)],
    influence_summary[, .(dgp, model_correct, n, median_max_rel_norm, q90_max_rel_norm, median_abs_top_delta_rho)],
    by = c("dgp", "model_correct", "n"),
    all = TRUE
  )
  setorder(manuscript, dgp, n)
  fwrite(manuscript, file.path(table_dir, "gof_influence_manuscript_summary.csv"))

  meta <- data.table(
    raw_dir = raw_dir,
    table_dir = table_dir,
    n_rows_raw = nrow(raw),
    n_files = length(list.files(raw_dir, pattern = "\\.rds$", full.names = FALSE)),
    B = unique(raw$B)[1L],
    generated_at = as.character(Sys.time())
  )
  fwrite(meta, file.path(table_dir, "gof_aggregation_metadata.csv"))

  log_message("GOF/influence aggregation finished.")
  log_message("Wrote tables to: ", table_dir)
  invisible(list(raw = raw, gof_summary = gof_summary, influence_summary = influence_summary))
}

if (isTRUE(aggregate_after)) {
  aggregate_gof_results(raw_dir = out_dir, table_dir = table_dir)
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
