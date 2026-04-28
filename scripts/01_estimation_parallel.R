#!/usr/bin/env Rscript
# scripts/01_estimation_parallel.R
# -----------------------------------------------------------------------------
# Parallel Monte Carlo experiment for point estimation and rho inference.
# Run from the project root, for example:
#   Rscript scripts/01_estimation_parallel.R --cores=14 --R=1500 --chunk=25
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

source("R/bs2_core.R")
source("R/sim_parallel_utils.R")
source("R/sim_designs.R")

cores <- get_int_arg("cores", max(1L, parallelly::availableCores() - 1L))
R <- get_int_arg("R", 1500L)
chunk_size <- get_int_arg("chunk", 25L)
seed <- get_int_arg("seed", 202501L)
backend <- get_arg("backend", "multicore")
design_name <- tolower(get_arg("design", "estimation"))
out_dir <- get_arg("out-dir", "results/raw_estimation")
prefix <- get_arg("prefix", paste0("est_", design_name))
force_cf <- get_arg("cf", NA_character_)

setup_parallel(cores = cores, backend = backend)
ensure_dirs(out_dir, "logs", "results/tables")

make_design <- function(name) {
  if (name %in% c("est", "estimation", "e1")) {
    make_estimation_design()
  } else if (name %in% c("rho", "rho_inference", "e2")) {
    make_rho_inference_design()
  } else if (name %in% c("debug", "test")) {
    make_debug_designs()$estimation
  } else {
    stop("Unknown --design value: ", name,
         ". Use estimation, rho_inference, or debug.")
  }
}

design <- make_design(design_name)

if (!is.na(force_cf)) {
  design[, do_cf := get_bool_arg("cf", TRUE)]
}
if (!"do_cf" %in% names(design)) {
  design[, do_cf := TRUE]
}
if (!"K_cf" %in% names(design)) {
  design[, K_cf := 5L]
}

tasks <- make_rep_tasks(
  design = design,
  R = R,
  chunk_size = chunk_size,
  out_dir = out_dir,
  prefix = prefix
)

tasks_todo <- tasks[!file.exists(out_file)]

message_header("BS2 estimation simulation")
message("Design: ", design_name)
message("Design cells: ", nrow(design))
message("Total tasks: ", nrow(tasks))
message("Pending tasks: ", nrow(tasks_todo))
message("R per cell: ", R)
message("Chunk size: ", chunk_size)
message("Workers: ", cores)
message("Backend: ", backend)
message("Output directory: ", out_dir)

one_estimation_replication <- function(cell, rep_id) {
  alpha <- c(cell$alpha1, cell$alpha2)
  beta <- c(cell$beta1, cell$beta2)
  rho <- if ("rho" %in% names(cell)) cell$rho else cell$rho_true
  n <- as.integer(cell$n)

  dat <- rbs2(n = n, alpha = alpha, beta = beta, rho = rho)

  fit_mm <- safe_eval(fit_bs2_mm(dat[, 1], dat[, 2]))

  fit_ml <- NULL
  if (fit_ok(fit_mm)) {
    fit_ml <- safe_eval(
      fit_bs2_profile(dat[, 1], dat[, 2], beta_start = fit_mm$beta)
    )
  }

  fit_pl <- NULL
  if (fit_ok(fit_mm)) {
    fit_pl <- safe_eval(
      fit_bs2_pl(
        dat[, 1], dat[, 2],
        alpha_hat = fit_mm$alpha,
        beta_hat = fit_mm$beta
      )
    )
  }

  fit_cf <- NULL
  if (isTRUE(cell$do_cf)) {
    fit_cf <- safe_eval(
      fit_bs2_pl_crossfit(
        dat[, 1], dat[, 2],
        K = as.integer(cell$K_cf),
        seed = NULL
      )
    )
  }

  data.table(
    rep = as.integer(rep_id),
    experiment = if ("experiment" %in% names(cell)) as.character(cell$experiment) else "E1_estimation",
    cell_id = if ("cell_id" %in% names(cell)) as.integer(cell$cell_id) else NA_integer_,
    alpha_id = if ("alpha_id" %in% names(cell)) as.character(cell$alpha_id) else NA_character_,
    n = n,
    alpha1_true = alpha[1],
    alpha2_true = alpha[2],
    beta1_true = beta[1],
    beta2_true = beta[2],
    rho_true = rho,

    mm_ok = fit_ok(fit_mm),
    ml_ok = fit_ok(fit_ml),
    pl_ok = fit_ok(fit_pl),
    cf_ok = fit_ok(fit_cf),

    mm_alpha1 = num_or_na(fit_mm, "alpha", 1),
    mm_alpha2 = num_or_na(fit_mm, "alpha", 2),
    mm_beta1  = num_or_na(fit_mm, "beta", 1),
    mm_beta2  = num_or_na(fit_mm, "beta", 2),
    mm_rho    = num_or_na(fit_mm, "rho"),

    ml_alpha1 = num_or_na(fit_ml, "alpha", 1),
    ml_alpha2 = num_or_na(fit_ml, "alpha", 2),
    ml_beta1  = num_or_na(fit_ml, "beta", 1),
    ml_beta2  = num_or_na(fit_ml, "beta", 2),
    ml_rho    = num_or_na(fit_ml, "rho"),
    ml_logLik = num_or_na(fit_ml, "logLik"),
    ml_convergence = num_or_na(fit_ml, "convergence"),

    pl_rho       = num_or_na(fit_pl, "rho"),
    pl_lambda    = num_or_na(fit_pl, "lambda"),
    pl_se_rho    = num_or_na(fit_pl, "se_rho"),
    pl_se_lambda = num_or_na(fit_pl, "se_lambda"),
    pl_lr0       = num_or_na(fit_pl, "lr0"),
    pl_score0    = num_or_na(fit_pl, "score0"),
    pl_separation = bool_or_na(fit_pl, "separation"),
    pl_boundary   = bool_or_na(fit_pl, "boundary"),

    cf_rho       = num_or_na(fit_cf, "rho"),
    cf_lambda    = num_or_na(fit_cf, "lambda"),
    cf_se_rho    = num_or_na(fit_cf, "se_rho"),
    cf_se_lambda = num_or_na(fit_cf, "se_lambda"),
    cf_lr0       = num_or_na(fit_cf, "lr0"),
    cf_score0    = num_or_na(fit_cf, "score0"),
    cf_separation = bool_or_na(fit_cf, "separation"),
    cf_boundary   = bool_or_na(fit_cf, "boundary")
  )
}

run_task <- function(task_row) {
  if (file.exists(task_row$out_file)) return(task_row$out_file)

  cell <- design[task_row$design_row]
  set.seed(task_seed(seed, task_row$cell_id, task_row$chunk_id))

  reps <- seq.int(task_row$rep_start, task_row$rep_end)
  ans <- rbindlist(
    lapply(reps, function(r) one_estimation_replication(cell, r)),
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
message("Estimation simulation finished.")
