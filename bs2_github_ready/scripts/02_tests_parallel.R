#!/usr/bin/env Rscript
# scripts/02_tests_parallel.R
# -----------------------------------------------------------------------------
# Parallel Monte Carlo experiment for empirical size and power of BS2 tests.
# Run from the project root, for example:
#   Rscript scripts/02_tests_parallel.R --design=size --cores=14 --R=3000
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

source("R/bs2_core.R")
source("R/sim_parallel_utils.R")
source("R/sim_designs.R")

cores <- get_int_arg("cores", max(1L, parallelly::availableCores() - 1L))
R <- get_int_arg("R", 3000L)
chunk_size <- get_int_arg("chunk", 25L)
seed <- get_int_arg("seed", 202502L)
backend <- get_arg("backend", "multicore")
design_name <- tolower(get_arg("design", "size"))
out_dir <- get_arg("out-dir", "results/raw_tests")
prefix <- get_arg("prefix", paste0("tests_", design_name))
level <- get_num_arg("level", 0.05)
B_sign <- get_int_arg("B-sign", 9999L)
sign_chunk <- get_int_arg("sign-chunk", 2000L)

setup_parallel(cores = cores, backend = backend)
ensure_dirs(out_dir, "logs", "results/tables")

make_design <- function(name) {
  if (name %in% c("size", "e3a")) {
    make_size_design()
  } else if (name %in% c("power", "e3b")) {
    make_power_design()
  } else if (name %in% c("debug", "test")) {
    make_debug_designs()$size
  } else {
    stop("Unknown --design value: ", name,
         ". Use size, power, or debug.")
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

message_header("BS2 test simulation")
message("Design: ", design_name)
message("Design cells: ", nrow(design))
message("Total tasks: ", nrow(tasks))
message("Pending tasks: ", nrow(tasks_todo))
message("R per cell: ", R)
message("Chunk size: ", chunk_size)
message("Workers: ", cores)
message("Backend: ", backend)
message("Level: ", level)
message("B-sign for H04 PL-score calibration: ", B_sign)
message("Output directory: ", out_dir)

one_test_replication <- function(cell, rep_id) {
  h <- as.character(cell$hypothesis)
  n <- as.integer(cell$n)
  alpha <- c(cell$alpha1, cell$alpha2)
  beta <- c(cell$beta1, cell$beta2)
  rho_true <- if ("rho_true" %in% names(cell)) cell$rho_true else cell$rho_design
  rho0 <- if ("rho0" %in% names(cell) && is.finite(cell$rho0)) cell$rho0 else 0.0
  df <- if ("df" %in% names(cell)) as.integer(cell$df) else restriction_rank(h)

  dat <- rbs2(n = n, alpha = alpha, beta = beta, rho = rho_true)

  fit_mm <- safe_eval(fit_bs2_mm(dat[, 1], dat[, 2]))

  fit_u <- NULL
  if (fit_ok(fit_mm)) {
    fit_u <- safe_eval(
      fit_bs2_profile(dat[, 1], dat[, 2], beta_start = fit_mm$beta)
    )
  }

  fit_r <- NULL
  if (fit_ok(fit_u)) {
    fit_r <- safe_eval(
      fit_bs2_restricted(
        dat[, 1], dat[, 2],
        hypothesis = h,
        rho0 = rho0,
        beta_start = fit_u$beta
      )
    )
  }

  tests <- NULL
  if (fit_ok(fit_u) && fit_ok(fit_r) && fit_ok(fit_mm)) {
    tests <- safe_eval(
      compute_bs2_tests(
        dat[, 1], dat[, 2],
        fit_u = fit_u,
        fit_r = fit_r,
        fit_mm = fit_mm,
        hypothesis = h,
        rho0 = rho0
      )
    )
  }

  crit <- qchisq(1 - level, df = df)
  lr <- num_or_na(tests, "LR")
  score <- num_or_na(tests, "Score")
  wald <- num_or_na(tests, "Wald")
  mm_wald <- num_or_na(tests, "MM_Wald")

  pl_lr <- NA_real_
  pl_score <- NA_real_
  pl_p_asym_lr <- NA_real_
  pl_p_asym_score <- NA_real_
  pl_sign_p <- NA_real_
  pl_ok <- FALSE

  if (h == "H04" && fit_ok(fit_mm)) {
    fit_pl <- safe_eval(
      fit_bs2_pl(
        dat[, 1], dat[, 2],
        alpha_hat = fit_mm$alpha,
        beta_hat = fit_mm$beta
      )
    )

    if (fit_ok(fit_pl)) {
      pl_ok <- TRUE
      pl_lr <- num_or_na(fit_pl, "lr0")
      pl_score <- num_or_na(fit_pl, "score0")
      pl_p_asym_lr <- pchisq(pl_lr, df = 1, lower.tail = FALSE)
      pl_p_asym_score <- pchisq(pl_score, df = 1, lower.tail = FALSE)

      if (B_sign > 0L) {
        tmp <- safe_eval(
          calibrate_pl_score_fast(
            x = fit_pl$x,
            c = fit_pl$c,
            B = B_sign,
            chunk = sign_chunk
          )
        )
        if (!inherits(tmp, "try-error")) pl_sign_p <- as.numeric(tmp)
      }
    }
  }

  data.table(
    rep = as.integer(rep_id),
    experiment = if ("experiment" %in% names(cell)) as.character(cell$experiment) else NA_character_,
    cell_id = if ("cell_id" %in% names(cell)) as.integer(cell$cell_id) else NA_integer_,
    hypothesis = h,
    alt_id = if ("alt_id" %in% names(cell)) as.character(cell$alt_id) else NA_character_,
    n = n,
    df = df,
    level = level,
    alpha1_true = alpha[1],
    alpha2_true = alpha[2],
    beta1_true = beta[1],
    beta2_true = beta[2],
    rho_true = rho_true,
    rho0 = rho0,

    fit_u_ok = fit_ok(fit_u),
    fit_r_ok = fit_ok(fit_r),
    fit_mm_ok = fit_ok(fit_mm),
    tests_ok = fit_ok(tests),
    pl_ok = pl_ok,

    LR = lr,
    Score = score,
    Wald = wald,
    MM_Wald = mm_wald,
    LR_p = num_or_na(tests, "LR_p"),
    Score_p = num_or_na(tests, "Score_p"),
    Wald_p = num_or_na(tests, "Wald_p"),
    MM_Wald_p = num_or_na(tests, "MM_Wald_p"),

    PL_LR = pl_lr,
    PL_Score = pl_score,
    PL_LR_p_asym = pl_p_asym_lr,
    PL_Score_p_asym = pl_p_asym_score,
    PL_sign_p = pl_sign_p,

    LR_reject = as.integer(!is.na(lr) && lr > crit),
    Score_reject = as.integer(!is.na(score) && score > crit),
    Wald_reject = as.integer(!is.na(wald) && wald > crit),
    MM_Wald_reject = as.integer(!is.na(mm_wald) && mm_wald > crit),
    PL_LR_reject = as.integer(!is.na(pl_lr) && pl_lr > qchisq(1 - level, df = 1)),
    PL_Score_reject = as.integer(!is.na(pl_score) && pl_score > qchisq(1 - level, df = 1)),
    PL_sign_reject = as.integer(!is.na(pl_sign_p) && pl_sign_p < level)
  )
}

run_task <- function(task_row) {
  if (file.exists(task_row$out_file)) return(task_row$out_file)

  cell <- design[task_row$design_row]
  set.seed(task_seed(seed, task_row$cell_id, task_row$chunk_id))

  reps <- seq.int(task_row$rep_start, task_row$rep_end)
  ans <- rbindlist(
    lapply(reps, function(r) one_test_replication(cell, r)),
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
message("Test simulation finished.")
