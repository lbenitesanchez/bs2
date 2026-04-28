#!/usr/bin/env Rscript
# scripts/03_bootstrap_lr_parallel.R
# -----------------------------------------------------------------------------
# Parallel bootstrap stress experiment for LR calibration under H01--H04.
# The outer Monte Carlo replications are parallelized. The inner bootstrap loop
# is sequential inside each worker to avoid nested parallelism.
# Run from project root, for example:
#   Rscript scripts/03_bootstrap_lr_parallel.R --cores=14 --R=500 --B=499
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

source("R/bs2_core.R")
source("R/sim_parallel_utils.R")
source("R/sim_designs.R")

cores <- get_int_arg("cores", max(1L, parallelly::availableCores() - 1L))
R <- get_int_arg("R", 500L)
B <- get_int_arg("B", 499L)
chunk_size <- get_int_arg("chunk", 5L)
seed <- get_int_arg("seed", 202503L)
backend <- get_arg("backend", "multicore")
design_name <- tolower(get_arg("design", "stress"))
out_dir <- get_arg("out-dir", "results/raw_bootstrap")
prefix <- get_arg("prefix", paste0("boot_", design_name))
level <- get_num_arg("level", 0.05)

setup_parallel(cores = cores, backend = backend)
ensure_dirs(out_dir, "logs", "results/tables")

make_design <- function(name) {
  if (name %in% c("stress", "bootstrap", "e4")) {
    make_bootstrap_stress_design()
  } else if (name %in% c("debug", "test")) {
    make_debug_designs()$bootstrap
  } else {
    stop("Unknown --design value: ", name,
         ". Use stress or debug.")
  }
}

design <- make_design(design_name)

tasks <- make_rep_tasks(
  design = design,
  R = R,
  chunk_size = chunk_size,
  out_dir = out_dir,
  prefix = prefix
)

tasks_todo <- tasks[!file.exists(out_file)]

message_header("BS2 bootstrap LR stress simulation")
message("Design: ", design_name)
message("Design cells: ", nrow(design))
message("Total tasks: ", nrow(tasks))
message("Pending tasks: ", nrow(tasks_todo))
message("Outer R per cell: ", R)
message("Inner B per outer replication: ", B)
message("Chunk size: ", chunk_size)
message("Workers: ", cores)
message("Backend: ", backend)
message("Output directory: ", out_dir)

one_bootstrap_replication <- function(cell, rep_id) {
  h <- as.character(cell$hypothesis)
  n <- as.integer(cell$n)
  alpha <- c(cell$alpha1, cell$alpha2)
  beta <- c(cell$beta1, cell$beta2)
  rho_true <- cell$rho_true
  rho0 <- if ("rho0" %in% names(cell) && is.finite(cell$rho0)) cell$rho0 else 0.0
  df <- if ("df" %in% names(cell)) as.integer(cell$df) else restriction_rank(h)

  dat <- rbs2(n = n, alpha = alpha, beta = beta, rho = rho_true)

  boot <- safe_eval(
    boot_lr_h0(
      dat[, 1], dat[, 2],
      hypothesis = h,
      B = B,
      rho0 = rho0,
      seed = NULL,
      parallel = FALSE,
      cores = 1L,
      verbose = FALSE
    )
  )

  lr_obs <- num_or_na(boot, "lr_obs")
  p_boot <- num_or_na(boot, "p_boot")
  p_asym <- num_or_na(boot, "p_asym")
  crit95 <- num_or_na(boot, "crit95")
  B_eff <- num_or_na(boot, "B_eff")
  failure_rate <- num_or_na(boot, "failure_rate")

  data.table(
    rep = as.integer(rep_id),
    experiment = if ("experiment" %in% names(cell)) as.character(cell$experiment) else "E4_bootstrap_LR",
    cell_id = if ("cell_id" %in% names(cell)) as.integer(cell$cell_id) else NA_integer_,
    scenario = if ("scenario" %in% names(cell)) as.character(cell$scenario) else NA_character_,
    hypothesis = h,
    n = n,
    df = df,
    level = level,
    alpha1_true = alpha[1],
    alpha2_true = alpha[2],
    beta1_true = beta[1],
    beta2_true = beta[2],
    rho_true = rho_true,
    rho0 = rho0,
    B = B,
    B_eff = B_eff,
    bootstrap_ok = fit_ok(boot),
    bootstrap_failure_rate = failure_rate,
    LR_obs = lr_obs,
    p_asym = p_asym,
    p_boot = p_boot,
    boot_crit95 = crit95,
    LR_asym_reject = as.integer(!is.na(p_asym) && p_asym < level),
    LR_boot_reject = as.integer(!is.na(p_boot) && p_boot < level)
  )
}

run_task <- function(task_row) {
  if (file.exists(task_row$out_file)) return(task_row$out_file)

  cell <- design[task_row$design_row]
  set.seed(task_seed(seed, task_row$cell_id, task_row$chunk_id))

  reps <- seq.int(task_row$rep_start, task_row$rep_end)
  ans <- rbindlist(
    lapply(reps, function(r) one_bootstrap_replication(cell, r)),
    fill = TRUE,
    use.names = TRUE
  )

  saveRDS(ans, task_row$out_file)
  task_row$out_file
}

if (nrow(tasks_todo) > 0L) {
  invisible(
    future_lapply(
      seq_len(nrow(tasks_todo)),
      function(ii) run_task(tasks_todo[ii]),
      future.seed = TRUE
    )
  )
}

reset_parallel()
message("Bootstrap LR stress simulation finished.")
