#!/usr/bin/env Rscript
# scripts/05_power_parallel.R
# -----------------------------------------------------------------------------
# Expanded BS2 empirical-power simulation for H01--H04.
#
# This script is the power analogue of scripts/02_tests_parallel.R.  The key
# difference is that data are generated from the *alternative* parameters stored
# in make_power_design(), while the restricted fit is still computed under the
# null hypothesis indicated by each row.
#
# Output:
#   - Raw checkpointed RDS files in results/raw_power_score_fixed/ by default.
#   - Aggregated CSV tables in results/tables_power_score_fixed/ by default.
#
# Typical formal run:
#   nohup Rscript scripts/05_power_parallel.R \
#     --R=1500 --chunk=25 --cores=${CORES} --backend=multicore \
#     --Bsign=9999 --seed=202504 \
#     --out=results/raw_power_score_fixed \
#     --table-dir=results/tables_power_score_fixed \
#     > logs/power_score_fixed.log 2>&1 &
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

source("R/bs2_core.R")
# The analytic projected-score implementation is integrated in R/bs2_core.R.
source("R/sim_parallel_utils.R")
source("R/sim_designs.R")
# -----------------------------------------------------------------------------
# Compatibility fallbacks
# -----------------------------------------------------------------------------
# Some earlier versions of R/sim_parallel_utils.R did not define these helpers.
# Define them here so this script works with both old and updated project folders.
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
    if (length(files) == 0L) return(data.table::data.table())
    data.table::rbindlist(lapply(files, readRDS), fill = TRUE)
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

# -----------------------------------------------------------------------------
# Command-line arguments
# -----------------------------------------------------------------------------

cores <- get_int_arg("cores", max(1L, parallelly::availableCores() - 1L))
R <- get_int_arg("R", 1500L)
chunk_size <- get_int_arg("chunk", 25L)
seed <- get_int_arg("seed", 202504L)
backend <- get_arg("backend", "multicore")
design_name <- tolower(get_arg("design", "power"))
debug <- get_bool_arg("debug", FALSE) || identical(design_name, "debug")
include_null <- get_bool_arg("include-null", TRUE)
aggregate_after <- get_bool_arg("aggregate", TRUE)

# Accept both --Bsign= and --B-sign= for convenience.
Bsign <- suppressWarnings(as.integer(get_arg("Bsign", get_arg("B-sign", "0"))))
if (is.na(Bsign)) stop("Argument --Bsign/--B-sign must be an integer.")

out_arg <- get_arg("out", get_arg("out-dir", NULL))
out_dir <- if (!is.null(out_arg)) {
  out_arg
} else if (debug) {
  "results/raw_power_debug"
} else {
  "results/raw_power_score_fixed"
}

table_arg <- get_arg("table-dir", NULL)
table_dir <- if (!is.null(table_arg)) {
  table_arg
} else if (debug) {
  "results/tables_power_debug"
} else {
  "results/tables_power_score_fixed"
}

prefix <- get_arg("prefix", if (debug) "power_debug" else "power")

setup_parallel(cores = cores, backend = backend)

# -----------------------------------------------------------------------------
# Design construction
# -----------------------------------------------------------------------------

if (debug) {
  design <- make_power_design(include_null = FALSE, n_grid_common = c(30L))
  design <- design[hypothesis %in% c("H01", "H04")]
  design <- design[seq_len(min(.N, 4L))]
  R <- min(R, 8L)
  chunk_size <- min(chunk_size, 4L)
  if (Bsign > 0L) Bsign <- min(Bsign, 99L)
} else {
  if (design_name %in% c("power", "e3b", "e3b_power")) {
    design <- make_power_design(include_null = include_null)
  } else if (design_name %in% c("power_alt", "alternatives", "alt")) {
    design <- make_power_design(include_null = FALSE)
  } else {
    stop("Unknown --design value: ", design_name,
         ". Use power, power_alt, or debug.")
  }
}

# Ensure cell_id is sequential after any debug subsetting.
design <- add_cell_id(design)

tasks <- make_rep_tasks(
  design = design,
  R = R,
  chunk_size = chunk_size,
  out_dir = out_dir,
  prefix = prefix
)

tasks_todo <- tasks[!file.exists(out_file)]

log_message("Power design cells: ", nrow(design))
log_message("Total tasks: ", nrow(tasks))
log_message("Pending tasks: ", nrow(tasks_todo))
log_message("R per cell: ", R)
log_message("Chunk size: ", chunk_size)
log_message("Workers: ", cores)
log_message("Backend: ", backend)
log_message("Bsign: ", Bsign)
log_message("Include null cases: ", include_null)
log_message("Output directory: ", out_dir)
log_message("Table directory: ", table_dir)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

reject_stat <- function(stat, crit) {
  if (is.finite(stat)) as.integer(stat > crit) else NA_integer_
}

reject_pvalue <- function(p, alpha = 0.05) {
  if (is.finite(p)) as.integer(p < alpha) else NA_integer_
}

is_null_case <- function(hypothesis, alpha, beta, rho_true, rho0) {
  hypothesis <- as.character(hypothesis)
  if (hypothesis == "H01") return(abs(alpha[1] - alpha[2]) < 1e-12)
  if (hypothesis == "H02") return(abs(beta[1] - beta[2]) < 1e-12)
  if (hypothesis == "H03") {
    return(abs(alpha[1] - alpha[2]) < 1e-12 && abs(beta[1] - beta[2]) < 1e-12)
  }
  if (hypothesis == "H04") return(abs(rho_true - rho0) < 1e-12)
  NA
}

one_power_replication <- function(cell, rep_id) {
  pin_numeric_threads()

  hypothesis <- as.character(cell$hypothesis)
  n <- as.integer(cell$n)
  df <- restriction_rank(hypothesis)

  alpha <- c(as.numeric(cell$alpha1), as.numeric(cell$alpha2))
  beta <- c(as.numeric(cell$beta1), as.numeric(cell$beta2))
  rho_true <- as.numeric(cell$rho_true)
  rho0 <- if ("rho0" %in% names(cell) && is.finite(cell$rho0)) as.numeric(cell$rho0) else 0.0
  true_null <- is_null_case(hypothesis, alpha = alpha, beta = beta,
                            rho_true = rho_true, rho0 = rho0)

  dat <- rbs2(n = n, alpha = alpha, beta = beta, rho = rho_true)
  t1 <- dat[, 1]
  t2 <- dat[, 2]

  fit_mm <- safe_eval(fit_bs2_mm(t1, t2))

  fit_u <- NULL
  if (fit_ok(fit_mm)) {
    fit_u <- safe_eval(fit_bs2_profile(t1, t2, beta_start = fit_mm$beta))
  }

  fit_r <- NULL
  if (fit_ok(fit_u)) {
    fit_r <- safe_eval(
      fit_bs2_restricted(
        t1, t2,
        hypothesis = hypothesis,
        rho0 = rho0,
        beta_start = fit_u$beta
      )
    )
  }

  tests <- NULL
  if (fit_ok(fit_u) && fit_ok(fit_r) && fit_ok(fit_mm)) {
    tests <- safe_eval(
      compute_bs2_tests(
        t1, t2,
        fit_u = fit_u,
        fit_r = fit_r,
        fit_mm = fit_mm,
        hypothesis = hypothesis,
        rho0 = rho0
      )
    )
  }

  crit <- qchisq(0.95, df = df)

  LR <- num_or_na(tests, "LR")
  Score <- num_or_na(tests, "Score")
  Wald <- num_or_na(tests, "Wald")
  MM_Wald <- num_or_na(tests, "MM_Wald")

  PL_LR <- NA_real_
  PL_Score <- NA_real_
  PL_sign_p <- NA_real_
  PL_ok <- NA

  # Conditional pseudo-likelihood tests only target H04: rho = rho0.
  # The formulas implemented in fit_bs2_pl() and calibrate_pl_score_fast()
  # are intended for the independence case rho0 = 0.
  if (hypothesis == "H04" && abs(rho0) < 1e-12 && fit_ok(fit_mm)) {
    fit_pl <- safe_eval(
      fit_bs2_pl(t1, t2, alpha_hat = fit_mm$alpha, beta_hat = fit_mm$beta)
    )
    PL_ok <- fit_ok(fit_pl)

    if (isTRUE(PL_ok)) {
      PL_LR <- num_or_na(fit_pl, "lr0")
      PL_Score <- num_or_na(fit_pl, "score0")

      if (Bsign > 0L) {
        tmp_p <- safe_eval(calibrate_pl_score_fast(x = fit_pl$x, c = fit_pl$c, B = Bsign))
        PL_sign_p <- if (fit_ok(tmp_p)) as.numeric(tmp_p) else NA_real_
      }
    }
  }

  data.table(
    rep = rep_id,
    experiment = ifelse("experiment" %in% names(cell), cell$experiment, "E3b_power"),
    cell_id = ifelse("cell_id" %in% names(cell), cell$cell_id, NA_integer_),
    hypothesis = hypothesis,
    alt_id = ifelse("alt_id" %in% names(cell), as.character(cell$alt_id), NA_character_),
    df = df,
    n = n,

    alpha1_true = alpha[1],
    alpha2_true = alpha[2],
    beta1_true = beta[1],
    beta2_true = beta[2],
    rho_true = rho_true,
    rho0 = ifelse(hypothesis == "H04", rho0, NA_real_),
    true_null = as.logical(true_null),

    alpha_ratio = ifelse("alpha_ratio" %in% names(cell), as.numeric(cell$alpha_ratio), alpha[2] / alpha[1]),
    beta_ratio = ifelse("beta_ratio" %in% names(cell), as.numeric(cell$beta_ratio), beta[2] / beta[1]),

    mm_ok = fit_ok(fit_mm),
    ml_ok = fit_ok(fit_u),
    restricted_ok = fit_ok(fit_r),
    tests_ok = fit_ok(tests),
    pl_ok = PL_ok,

    LR = LR,
    Score = Score,
    Wald = Wald,
    MM_Wald = MM_Wald,
    PL_LR = PL_LR,
    PL_Score = PL_Score,
    PL_sign_p = PL_sign_p,

    LR_reject = reject_stat(LR, crit),
    Score_reject = reject_stat(Score, crit),
    Wald_reject = reject_stat(Wald, crit),
    MM_Wald_reject = reject_stat(MM_Wald, crit),
    PL_LR_reject = reject_stat(PL_LR, qchisq(0.95, 1)),
    PL_Score_reject = reject_stat(PL_Score, qchisq(0.95, 1)),
    PL_sign_reject = reject_pvalue(PL_sign_p, alpha = 0.05)
  )
}

run_power_task <- function(task_row) {
  if (file.exists(task_row$out_file)) return(task_row$out_file)

  cell <- design[task_row$cell_id]
  set.seed(task_seed(seed, task_row$cell_id, task_row$chunk_id))
  reps <- seq.int(task_row$rep_start, task_row$rep_end)

  ans <- lapply(reps, function(r) one_power_replication(cell = cell, rep_id = r))
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
    function(ii) run_power_task(tasks_todo[ii]),
    future.seed = TRUE
  ))
}

log_message("Power simulation finished.")

# -----------------------------------------------------------------------------
# Aggregation
# -----------------------------------------------------------------------------

aggregate_power_results <- function(raw_dir, table_dir) {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  raw <- read_rds_dir(raw_dir)
  if (nrow(raw) == 0L) {
    warning("No RDS files found in ", raw_dir, ". No tables were written.")
    return(invisible(NULL))
  }

  reject_cols <- intersect(
    c("LR_reject", "Score_reject", "Wald_reject", "MM_Wald_reject",
      "PL_LR_reject", "PL_Score_reject", "PL_sign_reject"),
    names(raw)
  )

  group_cols <- intersect(
    c("experiment", "hypothesis", "alt_id", "df", "n",
      "alpha1_true", "alpha2_true", "beta1_true", "beta2_true",
      "rho_true", "rho0", "true_null", "alpha_ratio", "beta_ratio"),
    names(raw)
  )

  rejection <- rbindlist(lapply(reject_cols, function(rc) {
    x <- raw[!is.na(get(rc))]
    if (nrow(x) == 0L) return(NULL)
    stat_name <- sub("_reject$", "", rc)
    x[, .(
      R_eff = .N,
      rejection_rate = mean(get(rc)),
      mcse = sqrt(mean(get(rc)) * (1 - mean(get(rc))) / .N)
    ), by = group_cols][, statistic := stat_name][]
  }), fill = TRUE)

  if (nrow(rejection) > 0L) {
    setcolorder(rejection, c(group_cols, "statistic", "R_eff", "rejection_rate", "mcse"))
    setorder(rejection, hypothesis, alt_id, rho_true, n, statistic)
    fwrite(rejection, file.path(table_dir, "power_rejection_summary.csv"))

    # Convenient wide tables by hypothesis.
    for (h in sort(unique(rejection$hypothesis))) {
      tmp <- rejection[hypothesis == h]
      if (nrow(tmp) == 0L) next
      wide_cols <- intersect(
        c("hypothesis", "alt_id", "n", "rho_true", "rho0", "true_null",
          "alpha1_true", "alpha2_true", "beta1_true", "beta2_true",
          "alpha_ratio", "beta_ratio"),
        names(tmp)
      )
      wide <- dcast(
        tmp,
        as.formula(paste(paste(wide_cols, collapse = " + "), "~ statistic")),
        value.var = "rejection_rate"
      )
      setorder(wide, alt_id, rho_true, n)
      fwrite(wide, file.path(table_dir, paste0("power_pivot_", h, ".csv")))
    }
  }

  fit_cols <- intersect(c("mm_ok", "ml_ok", "restricted_ok", "tests_ok", "pl_ok"), names(raw))
  fit_rates <- rbindlist(lapply(fit_cols, function(fc) {
    x <- raw[!is.na(get(fc))]
    if (nrow(x) == 0L) return(NULL)
    x[, .(
      R_total = .N,
      success_rate = mean(as.logical(get(fc))),
      failure_rate = mean(!as.logical(get(fc)))
    ), by = group_cols][, component := fc][]
  }), fill = TRUE)

  if (nrow(fit_rates) > 0L) {
    setcolorder(fit_rates, c(group_cols, "component", "R_total", "success_rate", "failure_rate"))
    setorder(fit_rates, hypothesis, alt_id, rho_true, n, component)
    fwrite(fit_rates, file.path(table_dir, "power_fit_rates.csv"))
  }

  meta <- data.table(
    raw_dir = raw_dir,
    table_dir = table_dir,
    n_rows_raw = nrow(raw),
    n_files = length(list.files(raw_dir, pattern = "\\.rds$", full.names = FALSE)),
    generated_at = as.character(Sys.time())
  )
  fwrite(meta, file.path(table_dir, "power_aggregation_metadata.csv"))

  log_message("Power aggregation finished.")
  log_message("Wrote tables to: ", table_dir)
  invisible(list(raw = raw, rejection = rejection, fit_rates = fit_rates))
}

if (isTRUE(aggregate_after)) {
  aggregate_power_results(raw_dir = out_dir, table_dir = table_dir)
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
