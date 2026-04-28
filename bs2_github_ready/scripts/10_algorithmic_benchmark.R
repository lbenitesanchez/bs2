#!/usr/bin/env Rscript
# Algorithmic benchmark: profiled two-dimensional optimizer versus direct 5D optimizer.

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

source("R/bs2_core.R")
source("R/sim_parallel_utils.R")

cores <- get_int_arg("cores", max(1L, parallelly::availableCores() - 1L))
backend <- get_arg("backend", "multicore")
R <- get_int_arg("R", 200L)
seed <- get_int_arg("seed", 202507L)
out_dir <- get_arg("out-dir", "results/tables_algorithmic_benchmark")
raw_dir <- get_arg("raw-dir", "results/raw_algorithmic_benchmark")
chunk_size <- get_int_arg("chunk", 10L)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
setup_parallel(cores, backend)
on.exit(reset_parallel(), add = TRUE)

scenarios <- data.table(
  scenario = c("Moderate", "High-skew / strong dep."),
  n = c(250L, 250L),
  alpha1 = c(1.0, 1.5),
  alpha2 = c(1.5, 2.5),
  beta1 = c(1.0, 1.0),
  beta2 = c(1.5, 1.5),
  rho = c(0.60, 0.85)
)
scenarios[, cell_id := .I]

tasks <- make_rep_tasks(scenarios, R = R, chunk_size = chunk_size, out_dir = raw_dir, prefix = "benchmark")
tasks_todo <- tasks[!file.exists(out_file)]

score_eta_norm <- function(t1, t2, theta) {
  u <- score_bs2_analytic(t1, t2, theta)
  mult <- c(theta[1:4], 1 - theta[5]^2)
  sqrt(sum((u * mult)^2))
}

fit_profile_counted <- function(t1, t2, beta_start) {
  evals <- 0L
  loglik_prof <- function(eta_beta) {
    evals <<- evals + 1L
    beta <- exp(eta_beta)
    fit <- .bs2_profile_unrestricted_from_beta(t1, t2, beta = beta)
    if (is.null(fit)) return(-Inf)
    fit$logLik
  }
  objective <- function(eta_beta) {
    val <- loglik_prof(eta_beta)
    if (!is.finite(val)) return(.Machine$double.xmax / 100)
    -val
  }
  tm <- system.time({
    opt <- optim(log(beta_start), objective, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))
  })
  beta_hat <- exp(opt$par)
  fit <- .bs2_profile_unrestricted_from_beta(t1, t2, beta = beta_hat)
  theta <- c(fit$alpha, fit$beta, fit$rho)
  list(theta = theta, logLik = fit$logLik, convergence = opt$convergence,
       evals = evals, time_ms = unname(1000 * tm[["elapsed"]]), grad_norm = score_eta_norm(t1, t2, theta))
}

fit_full5d_counted <- function(t1, t2, start_theta) {
  evals <- 0L
  eta_start <- c(log(start_theta[1:4]), atanh(.clip_rho(start_theta[5], eps = 1e-8)))
  eta_to_theta <- function(eta) c(exp(eta[1:4]), tanh(eta[5]))
  objective <- function(eta) {
    evals <<- evals + 1L
    theta <- eta_to_theta(eta)
    val <- ell_bs2(t1, t2, alpha = theta[1:2], beta = theta[3:4], rho = theta[5])
    if (!is.finite(val)) return(.Machine$double.xmax / 100)
    -val
  }
  tm <- system.time({
    opt <- optim(eta_start, objective, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))
  })
  theta <- eta_to_theta(opt$par)
  list(theta = theta, logLik = -opt$value, convergence = opt$convergence,
       evals = evals, time_ms = unname(1000 * tm[["elapsed"]]), grad_norm = score_eta_norm(t1, t2, theta))
}

one_rep <- function(cell, rep_id) {
  dat <- rbs2(cell$n, alpha = c(cell$alpha1, cell$alpha2), beta = c(cell$beta1, cell$beta2), rho = cell$rho)
  mm <- fit_bs2_mm(dat[, 1], dat[, 2])
  theta_start <- c(mm$alpha, mm$beta, mm$rho)
  prof <- tryCatch(fit_profile_counted(dat[, 1], dat[, 2], beta_start = mm$beta), error = function(e) NULL)
  full <- tryCatch(fit_full5d_counted(dat[, 1], dat[, 2], start_theta = theta_start), error = function(e) NULL)
  rbindlist(list(
    data.table(rep_id = rep_id, scenario = cell$scenario, method = "Profile",
               ok = !is.null(prof), convergence = if (!is.null(prof)) prof$convergence else NA_integer_,
               logLik = if (!is.null(prof)) prof$logLik else NA_real_, evals = if (!is.null(prof)) prof$evals else NA_integer_,
               time_ms = if (!is.null(prof)) prof$time_ms else NA_real_, grad_norm = if (!is.null(prof)) prof$grad_norm else NA_real_),
    data.table(rep_id = rep_id, scenario = cell$scenario, method = "Full 5D",
               ok = !is.null(full), convergence = if (!is.null(full)) full$convergence else NA_integer_,
               logLik = if (!is.null(full)) full$logLik else NA_real_, evals = if (!is.null(full)) full$evals else NA_integer_,
               time_ms = if (!is.null(full)) full$time_ms else NA_real_, grad_norm = if (!is.null(full)) full$grad_norm else NA_real_)
  ), fill = TRUE)
}

run_task <- function(task_row) {
  if (file.exists(task_row$out_file)) return(task_row$out_file)
  cell <- scenarios[task_row$design_row]
  set.seed(task_seed(seed, task_row$cell_id, task_row$chunk_id))
  reps <- seq.int(task_row$rep_start, task_row$rep_end)
  ans <- rbindlist(lapply(reps, function(r) one_rep(cell, r)), fill = TRUE)
  saveRDS(ans, task_row$out_file)
  task_row$out_file
}

if (nrow(tasks_todo) > 0L) {
  invisible(future_lapply(seq_len(nrow(tasks_todo)), function(ii) run_task(tasks_todo[ii]), future.seed = TRUE))
}

raw <- read_rds_dir(raw_dir)
fwrite(raw, file.path(out_dir, "algorithmic_benchmark_raw.csv"))

summary <- raw[, .(
  stable_conv_pct = mean(ok & convergence == 0L & is.finite(grad_norm) & grad_norm <= 1e-4, na.rm = TRUE) * 100,
  median_evals = median(evals, na.rm = TRUE),
  median_time_ms = median(time_ms, na.rm = TRUE),
  median_grad_norm = median(grad_norm, na.rm = TRUE),
  median_logLik = median(logLik, na.rm = TRUE)
), by = .(scenario, method)]

# Add profile-vs-full log-likelihood gaps by scenario.
gaps <- dcast(raw[ok == TRUE, .(rep_id, scenario, method, logLik)], rep_id + scenario ~ method, value.var = "logLik")
if (all(c("Profile", "Full 5D") %in% names(gaps))) {
  gap_summary <- gaps[, .(
    median_abs_logLik_gap = median(abs(Profile - `Full 5D`), na.rm = TRUE),
    max_abs_logLik_gap = max(abs(Profile - `Full 5D`), na.rm = TRUE)
  ), by = scenario]
  summary <- merge(summary, gap_summary, by = "scenario", all.x = TRUE)
}

fwrite(summary, file.path(out_dir, "algorithmic_benchmark_summary.csv"))
cat("Algorithmic benchmark finished. Tables written to ", normalizePath(out_dir), "\n", sep = "")
