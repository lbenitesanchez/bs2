#!/usr/bin/env Rscript
# Reproduce the BMD empirical-application computations.
# Outputs CSV tables with unrestricted/restricted BS2 fits, tests, bootstrap LR
# p-values, pseudo-likelihood dependence summaries and a lognormal benchmark.

suppressPackageStartupMessages({
  library(data.table)
})

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(pat, "", hit[[1L]])
}
get_int_arg <- function(name, default) as.integer(get_arg(name, as.character(default)))
get_num_arg <- function(name, default) as.numeric(get_arg(name, as.character(default)))
get_bool_arg <- function(name, default = FALSE) {
  val <- tolower(as.character(get_arg(name, as.character(default))))
  val %in% c("true", "t", "yes", "y", "1")
}

bmd_file <- get_arg("data", "data/bmd.csv")
out_dir <- get_arg("out-dir", "results/tables")
B <- get_int_arg("B", 499L)
B_sign <- get_int_arg("B-sign", 99999L)
seed <- get_int_arg("seed", 202606L)
level <- get_num_arg("level", 0.05)
do_bootstrap <- get_bool_arg("bootstrap", TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("R/bs2_core.R")

if (!file.exists(bmd_file)) {
  stop("BMD data file not found: ", bmd_file,
       "\nRun from the repository root or pass --data=path/to/bmd.csv.", call. = FALSE)
}

bmd <- fread(bmd_file)
required <- c("id", "t1", "t2")
missing <- setdiff(required, names(bmd))
if (length(missing) > 0L) stop("Missing BMD columns: ", paste(missing, collapse = ", "))
bmd <- bmd[, .(id = as.integer(id), t1 = as.numeric(t1), t2 = as.numeric(t2))]

if (any(!is.finite(bmd$t1)) || any(!is.finite(bmd$t2)) || any(bmd$t1 <= 0) || any(bmd$t2 <= 0)) {
  stop("BMD variables must be finite and strictly positive.")
}

n <- nrow(bmd)
set.seed(seed)

# -----------------------------------------------------------------------------
# BS2 unrestricted ML/MM estimates and standard errors
# -----------------------------------------------------------------------------

fit_mm <- fit_bs2_mm(bmd$t1, bmd$t2)
fit_ml <- fit_bs2_profile(bmd$t1, bmd$t2, beta_start = fit_mm$beta)

I_ml <- fisher_info_bs2(fit_ml$alpha, fit_ml$beta, fit_ml$rho, gh_n = 25L)
Omega_ml <- .safe_solve(I_ml)
se_ml <- sqrt(diag(Omega_ml) / n)

Xi_mm <- mm_cov_bs2(fit_mm$alpha, fit_mm$beta, fit_mm$rho, gh_n = 25L)
se_mm <- sqrt(diag(Xi_mm) / n)

zcrit <- qnorm(1 - level / 2)
make_est_dt <- function(method, theta, se) {
  data.table(
    method = method,
    parameter = c("alpha1", "alpha2", "beta1", "beta2", "rho"),
    estimate = as.numeric(theta),
    se = as.numeric(se),
    ci_lower = as.numeric(theta) - zcrit * as.numeric(se),
    ci_upper = as.numeric(theta) + zcrit * as.numeric(se)
  )
}

estimates <- rbind(
  make_est_dt("ML", c(fit_ml$alpha, fit_ml$beta, fit_ml$rho), se_ml),
  make_est_dt("MM", c(fit_mm$alpha, fit_mm$beta, fit_mm$rho), se_mm)
)
fwrite(estimates, file.path(out_dir, "bmd_estimates.csv"))

# -----------------------------------------------------------------------------
# Restricted fits and H01--H04 tests
# -----------------------------------------------------------------------------

hypotheses <- c("H01", "H02", "H03", "H04")
restricted <- vector("list", length(hypotheses))
tests <- vector("list", length(hypotheses))

for (j in seq_along(hypotheses)) {
  h <- hypotheses[j]
  rho0 <- if (h == "H04") 0 else 0
  fit_r <- fit_bs2_restricted(bmd$t1, bmd$t2, hypothesis = h, rho0 = rho0, beta_start = fit_ml$beta)
  restricted[[j]] <- data.table(
    hypothesis = h,
    alpha1 = fit_r$alpha[1],
    alpha2 = fit_r$alpha[2],
    beta1 = fit_r$beta[1],
    beta2 = fit_r$beta[2],
    rho = fit_r$rho,
    logLik = fit_r$logLik,
    convergence = fit_r$convergence
  )

  tst <- compute_bs2_tests(bmd$t1, bmd$t2, fit_u = fit_ml, fit_r = fit_r, fit_mm = fit_mm, hypothesis = h, rho0 = rho0)
  boot_p <- NA_real_
  boot_crit <- NA_real_
  boot_fail <- NA_real_
  if (isTRUE(do_bootstrap) && B > 0L) {
    boot <- boot_lr_h0(bmd$t1, bmd$t2, hypothesis = h, B = B, rho0 = rho0, seed = seed + 100L + j,
                       beta_start = fit_ml$beta, parallel = FALSE)
    boot_p <- boot$p_boot
    boot_crit <- boot$crit95
    boot_fail <- boot$failure_rate
  }
  tests[[j]] <- data.table(
    hypothesis = h,
    df = tst$df,
    LR = tst$LR,
    LR_p = tst$LR_p,
    Score = tst$Score,
    Score_p = tst$Score_p,
    Wald = tst$Wald,
    Wald_p = tst$Wald_p,
    MM_Wald = tst$MM_Wald,
    MM_Wald_p = tst$MM_Wald_p,
    LR_boot_p = boot_p,
    LR_boot_crit95 = boot_crit,
    LR_boot_failure_rate = boot_fail
  )
}

fwrite(rbindlist(restricted), file.path(out_dir, "bmd_restricted_fits.csv"))
fwrite(rbindlist(tests), file.path(out_dir, "bmd_tests.csv"))

# -----------------------------------------------------------------------------
# Conditional pseudo-likelihood for dependence
# -----------------------------------------------------------------------------

fit_pl <- fit_bs2_pl(bmd$t1, bmd$t2, alpha_hat = fit_mm$alpha, beta_hat = fit_mm$beta)
set.seed(seed + 5000L)
pl_sign_p <- if (B_sign > 0L) calibrate_pl_score_fast(fit_pl$x, fit_pl$c, B = B_sign) else NA_real_

pl_interval <- function(fit, conf = 0.95, grid_n = 20001L) {
  rho_grid <- seq(-0.999, 0.999, length.out = grid_n)
  lambda_grid <- 2 * rho_grid / (1 - rho_grid^2)
  ell <- vapply(lambda_grid, function(lambda) {
    eta <- lambda * fit$x
    sum(fit$c * eta - log1pexp(eta))
  }, numeric(1L))
  keep <- 2 * (fit$logLik - ell) <= qchisq(conf, df = 1)
  if (!any(keep)) return(c(NA_real_, NA_real_))
  range(rho_grid[keep])
}
pl_ci <- pl_interval(fit_pl, conf = 1 - level)

pl_summary <- data.table(
  lambda = fit_pl$lambda,
  rho = fit_pl$rho,
  se_rho_working = fit_pl$se_rho,
  lr0 = fit_pl$lr0,
  lr0_p_chisq = pchisq(fit_pl$lr0, df = 1, lower.tail = FALSE),
  score0 = fit_pl$score0,
  score0_p_chisq = pchisq(fit_pl$score0, df = 1, lower.tail = FALSE),
  score0_p_sign = pl_sign_p,
  pl_ci_lower_working = pl_ci[1],
  pl_ci_upper_working = pl_ci[2],
  B_sign = B_sign
)
fwrite(pl_summary, file.path(out_dir, "bmd_pseudolikelihood.csv"))

# -----------------------------------------------------------------------------
# Bivariate lognormal benchmark
# -----------------------------------------------------------------------------

dmvn2_log <- function(y, mu, Sigma) {
  yc <- sweep(y, 2, mu, "-")
  invS <- solve(Sigma)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  q <- rowSums((yc %*% invS) * yc)
  -log(2 * pi) - 0.5 * logdet - 0.5 * q
}

fit_lognormal2 <- function(t1, t2) {
  y <- cbind(log(t1), log(t2))
  mu <- colMeans(y)
  yc <- sweep(y, 2, mu, "-")
  Sigma <- crossprod(yc) / nrow(y)
  logLik <- sum(dmvn2_log(y, mu, Sigma) - log(t1) - log(t2))
  list(mu = mu, Sigma = Sigma, logLik = logLik)
}

fit_ln <- fit_lognormal2(bmd$t1, bmd$t2)
model_comp <- data.table(
  model = c("BS2", "Bivariate lognormal"),
  logLik = c(fit_ml$logLik, fit_ln$logLik),
  k = c(5L, 5L)
)
model_comp[, `:=`(AIC = -2 * logLik + 2 * k, BIC = -2 * logLik + log(n) * k)]
fwrite(model_comp, file.path(out_dir, "bmd_model_comparison.csv"))

cat("BMD analysis finished. Tables written to: ", normalizePath(out_dir), "\n", sep = "")
cat("Unrestricted BS2 logLik: ", sprintf("%.6f", fit_ml$logLik), "\n", sep = "")
cat("Unrestricted BS2 theta: ", paste(sprintf("%.6f", c(fit_ml$alpha, fit_ml$beta, fit_ml$rho)), collapse = ", "), "\n", sep = "")
