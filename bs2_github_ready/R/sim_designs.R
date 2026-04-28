# R/sim_designs.R
# -----------------------------------------------------------------------------
# Design grids for the expanded Monte Carlo, testing, bootstrap and diagnostic
# experiments of the bivariate Birnbaum--Saunders (BS2) simulation study.
#
# This file contains only *design definitions*: parameter grids, null settings,
# alternative settings and small helper functions.  It intentionally contains no
# model fitting code.  The statistical routines should live in R/bs2_core.R and
# the parallel execution utilities in R/sim_parallel_utils.R.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it with install.packages('data.table').")
  }
  library(data.table)
})

# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

.valid_hypotheses <- c("H01", "H02", "H03", "H04")

validate_hypotheses <- function(hypothesis) {
  hypothesis <- as.character(hypothesis)
  bad <- setdiff(hypothesis, .valid_hypotheses)
  if (length(bad) > 0L) {
    stop("Unknown hypothesis label(s): ", paste(bad, collapse = ", "),
         ". Valid labels are: ", paste(.valid_hypotheses, collapse = ", "), ".")
  }
  hypothesis
}

restriction_rank <- function(hypothesis) {
  hypothesis <- validate_hypotheses(hypothesis)
  vapply(
    hypothesis,
    function(h) {
      switch(h,
             H01 = 1L,
             H02 = 1L,
             H03 = 2L,
             H04 = 1L)
    },
    integer(1L)
  )
}

hypothesis_description <- function(hypothesis) {
  hypothesis <- validate_hypotheses(hypothesis)
  vapply(
    hypothesis,
    function(h) {
      switch(h,
             H01 = "alpha1 = alpha2",
             H02 = "beta1 = beta2",
             H03 = "alpha1 = alpha2 and beta1 = beta2",
             H04 = "rho = rho0")
    },
    character(1L)
  )
}

as_theta_list <- function(alpha1, alpha2, beta1, beta2, rho) {
  list(
    alpha = c(as.numeric(alpha1), as.numeric(alpha2)),
    beta  = c(as.numeric(beta1), as.numeric(beta2)),
    rho   = as.numeric(rho)
  )
}

theta_from_design_row <- function(row, rho_col = "rho") {
  if (!is.data.frame(row) || nrow(row) != 1L) {
    stop("row must be a one-row data.frame or data.table.")
  }
  if (!rho_col %in% names(row)) {
    stop("Column '", rho_col, "' was not found in row.")
  }
  as_theta_list(
    alpha1 = row$alpha1,
    alpha2 = row$alpha2,
    beta1  = row$beta1,
    beta2  = row$beta2,
    rho    = row[[rho_col]]
  )
}

add_cell_id <- function(design) {
  design <- as.data.table(design)
  if ("cell_id" %in% names(design)) design[, cell_id := NULL]
  design[, cell_id := .I]
  setcolorder(design, c("cell_id", setdiff(names(design), "cell_id")))
  design[]
}

# -----------------------------------------------------------------------------
# Global parameter grids
# -----------------------------------------------------------------------------

make_alpha_grid <- function() {
  data.table(
    alpha_id = c("mild", "baseline", "severe"),
    alpha1   = c(0.3, 1.0, 1.5),
    alpha2   = c(0.3, 1.5, 2.5),
    alpha_description = c(
      "mild marginal skewness",
      "baseline configuration",
      "severe marginal skewness"
    )
  )
}

make_estimation_rho_grid <- function() {
  c(-0.75, -0.30, 0.00, 0.30, 0.60, 0.85)
}

make_rho_inference_grid <- function() {
  c(-0.85, -0.60, -0.30, 0.00, 0.30, 0.60, 0.85)
}

make_size_rho_grid <- function() {
  c(-0.60, 0.00, 0.60)
}

# -----------------------------------------------------------------------------
# E1. Estimation grid: MM vs ML over skewness and dependence levels
# -----------------------------------------------------------------------------

make_estimation_design <- function(
    alpha_ids = c("mild", "baseline", "severe"),
    rho_grid = make_estimation_rho_grid(),
    n_grid = c(30L, 60L, 120L, 250L),
    beta = c(1.0, 1.5),
    do_cf = TRUE) {

  alpha_grid <- make_alpha_grid()[alpha_id %in% alpha_ids]
  if (nrow(alpha_grid) == 0L) {
    stop("No valid alpha_id selected. Available alpha_id values are: ",
         paste(make_alpha_grid()$alpha_id, collapse = ", "), ".")
  }
  if (length(beta) != 2L || any(!is.finite(beta)) || any(beta <= 0)) {
    stop("beta must be a numeric vector of length 2 with positive entries.")
  }

  design <- CJ(
    alpha_id = alpha_grid$alpha_id,
    rho = as.numeric(rho_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design <- merge(design, alpha_grid, by = "alpha_id", sort = FALSE)
  design[, `:=`(
    beta1 = as.numeric(beta[1L]),
    beta2 = as.numeric(beta[2L]),
    do_cf = as.logical(do_cf),
    experiment = "E1_estimation"
  )]

  setorder(design, alpha_id, rho, n)
  setcolorder(design, c(
    "experiment", "alpha_id", "alpha_description", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho", "do_cf"
  ))

  add_cell_id(design)
}

# -----------------------------------------------------------------------------
# E2. Dependence-parameter inference grid: ML/MM/PL/CF for rho
# -----------------------------------------------------------------------------

make_rho_inference_design <- function(
    alpha_ids = c("mild", "baseline", "severe"),
    rho_grid = make_rho_inference_grid(),
    n_grid = c(30L, 60L, 120L),
    beta = c(1.0, 1.5),
    K_cf = 5L) {

  alpha_grid <- make_alpha_grid()[alpha_id %in% alpha_ids]
  if (nrow(alpha_grid) == 0L) {
    stop("No valid alpha_id selected. Available alpha_id values are: ",
         paste(make_alpha_grid()$alpha_id, collapse = ", "), ".")
  }
  if (length(beta) != 2L || any(!is.finite(beta)) || any(beta <= 0)) {
    stop("beta must be a numeric vector of length 2 with positive entries.")
  }

  design <- CJ(
    alpha_id = alpha_grid$alpha_id,
    rho = as.numeric(rho_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design <- merge(design, alpha_grid, by = "alpha_id", sort = FALSE)
  design[, `:=`(
    beta1 = as.numeric(beta[1L]),
    beta2 = as.numeric(beta[2L]),
    K_cf = as.integer(K_cf),
    experiment = "E2_rho_inference"
  )]

  setorder(design, alpha_id, rho, n)
  setcolorder(design, c(
    "experiment", "alpha_id", "alpha_description", "n", "K_cf",
    "alpha1", "alpha2", "beta1", "beta2", "rho"
  ))

  add_cell_id(design)
}

# -----------------------------------------------------------------------------
# Null parameter settings for H01--H04
# -----------------------------------------------------------------------------

make_null_parameters <- function(hypothesis, rho = NULL, rho_design = NULL, rho0 = 0) {
  hypothesis <- validate_hypotheses(hypothesis)
  if (length(hypothesis) != 1L) {
    stop("make_null_parameters() expects a single hypothesis label.")
  }

  if (is.null(rho_design)) {
    rho_design <- if (!is.null(rho)) rho else rho0
  }
  rho_design <- as.numeric(rho_design)
  rho0 <- as.numeric(rho0)

  if (hypothesis == "H01") {
    out <- list(alpha = c(1.0, 1.0), beta = c(1.0, 1.5), rho = rho_design)
  } else if (hypothesis == "H02") {
    out <- list(alpha = c(1.0, 1.5), beta = c(1.0, 1.0), rho = rho_design)
  } else if (hypothesis == "H03") {
    out <- list(alpha = c(1.0, 1.0), beta = c(1.0, 1.0), rho = rho_design)
  } else if (hypothesis == "H04") {
    out <- list(alpha = c(1.0, 1.5), beta = c(1.0, 1.5), rho = rho0)
  }

  out$hypothesis <- hypothesis
  out$rho0 <- if (hypothesis == "H04") rho0 else NA_real_
  out$df <- restriction_rank(hypothesis)
  out$description <- hypothesis_description(hypothesis)
  out
}

# -----------------------------------------------------------------------------
# E3a. Empirical size grid for H01--H04
# -----------------------------------------------------------------------------

make_size_design <- function(
    hypotheses = c("H01", "H02", "H03", "H04"),
    rho_grid = make_size_rho_grid(),
    n_grid = c(24L, 30L, 60L, 120L),
    h04_rho0_grid = 0.0) {

  hypotheses <- validate_hypotheses(hypotheses)

  parts <- vector("list", length(hypotheses))
  names(parts) <- hypotheses

  for (h in hypotheses) {
    if (h == "H04") {
      tmp <- CJ(
        hypothesis = h,
        rho_design = as.numeric(h04_rho0_grid),
        rho0 = as.numeric(h04_rho0_grid),
        n = as.integer(n_grid),
        sorted = FALSE
      )
    } else {
      tmp <- CJ(
        hypothesis = h,
        rho_design = as.numeric(rho_grid),
        rho0 = NA_real_,
        n = as.integer(n_grid),
        sorted = FALSE
      )
    }
    parts[[h]] <- tmp
  }

  design <- rbindlist(parts, use.names = TRUE, fill = TRUE)

  design[, `:=`(
    alpha1 = NA_real_, alpha2 = NA_real_,
    beta1  = NA_real_, beta2  = NA_real_,
    rho_true = NA_real_,
    df = restriction_rank(hypothesis),
    null_description = hypothesis_description(hypothesis),
    experiment = "E3a_size"
  )]

  for (i in seq_len(nrow(design))) {
    par0 <- make_null_parameters(
      hypothesis = design$hypothesis[i],
      rho_design = design$rho_design[i],
      rho0 = ifelse(is.na(design$rho0[i]), 0.0, design$rho0[i])
    )
    design$alpha1[i] <- par0$alpha[1L]
    design$alpha2[i] <- par0$alpha[2L]
    design$beta1[i]  <- par0$beta[1L]
    design$beta2[i]  <- par0$beta[2L]
    design$rho_true[i] <- par0$rho
  }

  setorder(design, hypothesis, rho_design, n)
  setcolorder(design, c(
    "experiment", "hypothesis", "null_description", "df", "n",
    "alpha1", "alpha2", "beta1", "beta2",
    "rho_true", "rho_design", "rho0"
  ))

  add_cell_id(unique(design))
}

# -----------------------------------------------------------------------------
# E3b. Power grids for H01, H02, H03 and H04
# -----------------------------------------------------------------------------

make_power_design_H01 <- function(
    alpha_ratio_grid = c(1.00, 1.25, 1.50, 2.00),
    rho_grid = c(0.0, 0.6),
    n_grid = c(30L, 60L, 120L),
    beta = c(1.0, 1.5)) {

  design <- CJ(
    alpha_ratio = as.numeric(alpha_ratio_grid),
    rho_true = as.numeric(rho_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design[, `:=`(
    experiment = "E3b_power",
    hypothesis = "H01",
    alt_id = sprintf("alpha_ratio_%s", format(alpha_ratio, trim = TRUE)),
    alpha1 = 1.0,
    alpha2 = alpha_ratio,
    beta1 = as.numeric(beta[1L]),
    beta2 = as.numeric(beta[2L]),
    rho0 = NA_real_,
    df = 1L
  )]

  setcolorder(design, c(
    "experiment", "hypothesis", "alt_id", "df", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho_true", "rho0", "alpha_ratio"
  ))

  design[]
}

make_power_design_H02 <- function(
    beta_ratio_grid = c(1.00, 1.25, 1.50, 2.00),
    rho_grid = c(0.0, 0.6),
    n_grid = c(30L, 60L, 120L),
    alpha = c(1.0, 1.5)) {

  design <- CJ(
    beta_ratio = as.numeric(beta_ratio_grid),
    rho_true = as.numeric(rho_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design[, `:=`(
    experiment = "E3b_power",
    hypothesis = "H02",
    alt_id = sprintf("beta_ratio_%s", format(beta_ratio, trim = TRUE)),
    alpha1 = as.numeric(alpha[1L]),
    alpha2 = as.numeric(alpha[2L]),
    beta1 = 1.0,
    beta2 = beta_ratio,
    rho0 = NA_real_,
    df = 1L
  )]

  setcolorder(design, c(
    "experiment", "hypothesis", "alt_id", "df", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho_true", "rho0", "beta_ratio"
  ))

  design[]
}

make_power_design_H03 <- function(
    rho_grid = c(0.0, 0.6),
    n_grid = c(30L, 60L, 120L)) {

  alternatives <- data.table(
    alt_id = c("shape_only", "scale_only", "shape_and_scale"),
    alpha1 = c(1.0, 1.0, 1.0),
    alpha2 = c(1.5, 1.0, 1.5),
    beta1  = c(1.0, 1.0, 1.0),
    beta2  = c(1.0, 1.5, 1.5)
  )

  design <- CJ(
    alt_id = alternatives$alt_id,
    rho_true = as.numeric(rho_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design <- merge(design, alternatives, by = "alt_id", sort = FALSE)
  design[, `:=`(
    experiment = "E3b_power",
    hypothesis = "H03",
    rho0 = NA_real_,
    df = 2L
  )]

  setorder(design, alt_id, rho_true, n)
  setcolorder(design, c(
    "experiment", "hypothesis", "alt_id", "df", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho_true", "rho0"
  ))

  design[]
}

make_power_design_H04 <- function(
    rho_alt_grid = c(-0.60, -0.30, 0.00, 0.30, 0.60),
    rho0 = 0.0,
    n_grid = c(24L, 30L, 60L, 120L),
    alpha = c(1.0, 1.5),
    beta = c(1.0, 1.5)) {

  design <- CJ(
    rho_true = as.numeric(rho_alt_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design[, `:=`(
    experiment = "E3b_power",
    hypothesis = "H04",
    alt_id = sprintf("rho_%s", format(rho_true, trim = TRUE)),
    alpha1 = as.numeric(alpha[1L]),
    alpha2 = as.numeric(alpha[2L]),
    beta1 = as.numeric(beta[1L]),
    beta2 = as.numeric(beta[2L]),
    rho0 = as.numeric(rho0),
    df = 1L
  )]

  setcolorder(design, c(
    "experiment", "hypothesis", "alt_id", "df", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho_true", "rho0"
  ))

  design[]
}

make_power_design <- function(
    include_null = TRUE,
    n_grid_common = c(30L, 60L, 120L)) {

  out <- rbindlist(
    list(
      make_power_design_H01(n_grid = n_grid_common),
      make_power_design_H02(n_grid = n_grid_common),
      make_power_design_H03(n_grid = n_grid_common),
      make_power_design_H04(n_grid = c(24L, n_grid_common))
    ),
    use.names = TRUE,
    fill = TRUE
  )

  if (!isTRUE(include_null)) {
    out <- out[!(hypothesis == "H01" & alpha1 == alpha2)]
    out <- out[!(hypothesis == "H02" & beta1 == beta2)]
    out <- out[!(hypothesis == "H04" & abs(rho_true - rho0) < 1e-12)]
  }

  setorder(out, hypothesis, alt_id, rho_true, n)
  add_cell_id(out)
}

# -----------------------------------------------------------------------------
# E4. Bootstrap LR stress grid
# -----------------------------------------------------------------------------

make_bootstrap_stress_design <- function() {
  design <- data.table(
    scenario = paste0("B", 1:6),
    hypothesis = c("H01", "H02", "H03", "H04", "H04", "H04"),
    n = c(24L, 24L, 24L, 24L, 30L, 60L),
    alpha1 = c(1.0, 1.0, 1.0, 1.5, 1.5, 1.5),
    alpha2 = c(1.0, 1.5, 1.0, 2.5, 2.5, 2.5),
    beta1  = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
    beta2  = c(1.5, 1.0, 1.0, 1.5, 1.5, 1.5),
    rho_true = c(0.6, 0.6, 0.6, 0.0, 0.0, 0.0),
    rho0 = c(NA_real_, NA_real_, NA_real_, 0.0, 0.0, 0.0)
  )

  design[, `:=`(
    experiment = "E4_bootstrap_LR",
    df = restriction_rank(hypothesis),
    null_description = hypothesis_description(hypothesis)
  )]

  setcolorder(design, c(
    "experiment", "scenario", "hypothesis", "null_description", "df", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho_true", "rho0"
  ))

  add_cell_id(design)
}

# -----------------------------------------------------------------------------
# E5. Goodness-of-fit and influence diagnostic design
# -----------------------------------------------------------------------------

make_gof_design <- function(
    n_grid = c(24L, 50L, 100L),
    dgp_grid = c("bs2", "lognormal", "contaminated_bs2"),
    alpha = c(0.1491, 0.1674),
    beta = c(0.8313, 0.8292),
    rho = 0.9343,
    eps_contam = 0.05,
    beta_multiplier_contam = c(0.6, 0.6)) {

  design <- CJ(
    dgp = as.character(dgp_grid),
    n = as.integer(n_grid),
    sorted = FALSE
  )

  design[, `:=`(
    experiment = "E5_gof_influence",
    alpha1 = as.numeric(alpha[1L]),
    alpha2 = as.numeric(alpha[2L]),
    beta1 = as.numeric(beta[1L]),
    beta2 = as.numeric(beta[2L]),
    rho_true = as.numeric(rho),
    eps_contam = as.numeric(eps_contam),
    beta1_contam = as.numeric(beta[1L] * beta_multiplier_contam[1L]),
    beta2_contam = as.numeric(beta[2L] * beta_multiplier_contam[2L])
  )]

  setorder(design, dgp, n)
  setcolorder(design, c(
    "experiment", "dgp", "n",
    "alpha1", "alpha2", "beta1", "beta2", "rho_true",
    "eps_contam", "beta1_contam", "beta2_contam"
  ))

  add_cell_id(design)
}

# -----------------------------------------------------------------------------
# BMD data used in the real-data illustration
# -----------------------------------------------------------------------------

make_bmd_data <- function() {
  data.table(
    id = seq_len(24L),
    t1 = c(
      1.103, 0.842, 0.925, 0.857, 0.795, 0.787,
      0.933, 0.799, 0.945, 0.921, 0.792, 0.815,
      0.755, 0.880, 0.900, 0.764, 0.733, 0.932,
      0.856, 0.890, 0.688, 0.940, 0.493, 0.835
    ),
    t2 = c(
      1.027, 0.857, 0.875, 0.873, 0.811, 0.640,
      0.947, 0.886, 0.991, 0.977, 0.825, 0.851,
      0.770, 0.912, 0.905, 0.756, 0.765, 0.932,
      0.843, 0.879, 0.673, 0.949, 0.463, 0.776
    )
  )
}

# -----------------------------------------------------------------------------
# Debug/minimal designs for quick smoke tests
# -----------------------------------------------------------------------------

make_debug_designs <- function() {
  list(
    estimation = make_estimation_design(
      alpha_ids = "baseline",
      rho_grid = c(0.0, 0.6),
      n_grid = c(30L),
      do_cf = TRUE
    ),
    rho_inference = make_rho_inference_design(
      alpha_ids = "baseline",
      rho_grid = c(0.0, 0.6),
      n_grid = c(30L),
      K_cf = 5L
    ),
    size = make_size_design(
      hypotheses = c("H01", "H04"),
      rho_grid = c(0.0),
      n_grid = c(24L, 30L),
      h04_rho0_grid = 0.0
    ),
    bootstrap = make_bootstrap_stress_design()[scenario %in% c("B1", "B4")],
    gof = make_gof_design(n_grid = c(24L), dgp_grid = c("bs2"))
  )
}

make_all_designs <- function() {
  list(
    E1_estimation = make_estimation_design(),
    E2_rho_inference = make_rho_inference_design(),
    E3a_size = make_size_design(),
    E3b_power = make_power_design(),
    E4_bootstrap_LR = make_bootstrap_stress_design(),
    E5_gof_influence = make_gof_design(),
    BMD_data = make_bmd_data()
  )
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
