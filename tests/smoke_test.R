#!/usr/bin/env Rscript
# Short numerical smoke test for the BS2 reproducibility bundle.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/bs2_core.R")

stopifnot(file.exists("data/bmd.csv"))
bmd <- fread("data/bmd.csv")
stopifnot(all(c("id", "t1", "t2") %in% names(bmd)))

mm <- fit_bs2_mm(bmd$t1, bmd$t2)
ml <- fit_bs2_profile(bmd$t1, bmd$t2, beta_start = mm$beta)
stopifnot(is.finite(ml$logLik), abs(ml$rho) < 1, all(ml$alpha > 0), all(ml$beta > 0))

# Tolerances are intentionally loose enough to allow small platform differences.
stopifnot(abs(ml$logLik - 54.244) < 0.05)
stopifnot(abs(ml$rho - 0.9343) < 0.02)

fit_r <- fit_bs2_restricted(bmd$t1, bmd$t2, hypothesis = "H04", rho0 = 0, beta_start = ml$beta)
tst <- compute_bs2_tests(bmd$t1, bmd$t2, fit_u = ml, fit_r = fit_r, fit_mm = mm, hypothesis = "H04", rho0 = 0)
stopifnot(is.finite(tst$LR), tst$LR > 0, is.finite(tst$Score), is.finite(tst$Wald))

set.seed(1)
x <- rbs2(12, alpha = c(1, 1.5), beta = c(1, 1.5), rho = 0.4)
mm2 <- fit_bs2_mm(x[, 1], x[, 2])
ml2 <- fit_bs2_profile(x[, 1], x[, 2], beta_start = mm2$beta)
pl2 <- fit_bs2_pl(x[, 1], x[, 2], alpha_hat = mm2$alpha, beta_hat = mm2$beta)
stopifnot(is.finite(ml2$logLik), all(ml2$alpha > 0), all(ml2$beta > 0), abs(ml2$rho) < 1)
stopifnot(is.finite(pl2$rho) || isTRUE(pl2$boundary))

cat("Smoke test passed.\n")
