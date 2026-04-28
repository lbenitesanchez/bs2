#!/usr/bin/env Rscript
# Supplementary BMD contour comparison: BS2 versus bivariate lognormal.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(pat, "", hit[[1L]])
}

bmd_file <- get_arg("data", "data/bmd.csv")
fig_dir <- get_arg("fig-dir", "figures")
seed <- as.integer(get_arg("seed", "20260427"))
nsim_hdr <- as.integer(get_arg("nsim-hdr", "200000"))
formats <- tolower(trimws(unlist(strsplit(get_arg("formats", "pdf,png"), ",", fixed = TRUE))))
formats <- formats[nzchar(formats)]

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
source("R/bs2_core.R")

bmd <- fread(bmd_file)
stopifnot(all(c("id", "t1", "t2") %in% names(bmd)))
bmd <- bmd[, .(id = as.integer(id), t1 = as.numeric(t1), t2 = as.numeric(t2))]

fit_mm <- fit_bs2_mm(bmd$t1, bmd$t2)
fit_bs2 <- fit_bs2_profile(bmd$t1, bmd$t2, beta_start = fit_mm$beta)

# Vectorized BS2 density for contour computation.
dbs2_mat <- function(x, alpha, beta, rho, log = FALSE) {
  x <- as.matrix(x)
  out <- rep(-Inf, nrow(x))
  ok <- x[, 1] > 0 & x[, 2] > 0 & all(alpha > 0) & all(beta > 0) & abs(rho) < 1
  if (any(ok)) {
    t1 <- x[ok, 1]
    t2 <- x[ok, 2]
    z1 <- a_bs(t1, alpha[1], beta[1])
    z2 <- a_bs(t2, alpha[2], beta[2])
    q <- (z1^2 - 2 * rho * z1 * z2 + z2^2) / (1 - rho^2)
    out[ok] <- -log(2 * pi) - 0.5 * log(1 - rho^2) - 0.5 * q +
      log(A_bs(t1, alpha[1], beta[1])) + log(A_bs(t2, alpha[2], beta[2]))
  }
  if (log) out else exp(out)
}

# Bivariate lognormal benchmark.
dmvn2_log <- function(y, mu, Sigma) {
  y <- as.matrix(y)
  yc <- sweep(y, 2, mu, "-")
  invS <- solve(Sigma)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  q <- rowSums((yc %*% invS) * yc)
  -log(2 * pi) - 0.5 * logdet - 0.5 * q
}

dlnorm2 <- function(x, mu, Sigma, log = FALSE) {
  x <- as.matrix(x)
  out <- rep(-Inf, nrow(x))
  ok <- x[, 1] > 0 & x[, 2] > 0
  if (any(ok)) {
    lx <- log(x[ok, , drop = FALSE])
    out[ok] <- dmvn2_log(lx, mu, Sigma) - log(x[ok, 1]) - log(x[ok, 2])
  }
  if (log) out else exp(out)
}

rlnorm2 <- function(n, mu, Sigma) {
  R <- chol(Sigma)
  z <- matrix(rnorm(2 * n), nrow = n, ncol = 2) %*% R
  exp(sweep(z, 2, mu, "+"))
}

fit_lognormal2 <- function(t1, t2) {
  y <- cbind(log(t1), log(t2))
  mu <- colMeans(y)
  yc <- sweep(y, 2, mu, "-")
  Sigma <- crossprod(yc) / nrow(y)
  list(mu = mu, Sigma = Sigma,
       logLik = sum(dlnorm2(cbind(t1, t2), mu = mu, Sigma = Sigma, log = TRUE)))
}

fit_ln <- fit_lognormal2(bmd$t1, bmd$t2)

set.seed(seed)
hdr_probs <- c(0.50, 0.75, 0.90, 0.95)
sim_bs2 <- rbs2(nsim_hdr, fit_bs2$alpha, fit_bs2$beta, fit_bs2$rho)
sim_ln <- rlnorm2(nsim_hdr, fit_ln$mu, fit_ln$Sigma)

hdr_thresholds <- function(xsim, dens_fun, probs) {
  dens <- dens_fun(xsim)
  as.numeric(quantile(dens, probs = 1 - probs, type = 8, names = FALSE))
}
thr_bs2 <- hdr_thresholds(sim_bs2, function(x) dbs2_mat(x, fit_bs2$alpha, fit_bs2$beta, fit_bs2$rho), hdr_probs)
thr_ln <- hdr_thresholds(sim_ln, function(x) dlnorm2(x, fit_ln$mu, fit_ln$Sigma), hdr_probs)

pad_range <- function(r, frac = 0.08) {
  w <- diff(r)
  c(max(1e-6, r[1] - frac * w), r[2] + frac * w)
}
x_range <- range(bmd$t1, quantile(sim_bs2[, 1], c(0.0025, 0.9975)), quantile(sim_ln[, 1], c(0.0025, 0.9975)))
y_range <- range(bmd$t2, quantile(sim_bs2[, 2], c(0.0025, 0.9975)), quantile(sim_ln[, 2], c(0.0025, 0.9975)))
xg <- seq(pad_range(x_range)[1], pad_range(x_range)[2], length.out = 300)
yg <- seq(pad_range(y_range)[1], pad_range(y_range)[2], length.out = 300)
grid <- expand.grid(t1 = xg, t2 = yg)

z_bs2 <- matrix(dbs2_mat(as.matrix(grid), fit_bs2$alpha, fit_bs2$beta, fit_bs2$rho), nrow = length(xg), ncol = length(yg))
z_ln <- matrix(dlnorm2(as.matrix(grid), fit_ln$mu, fit_ln$Sigma), nrow = length(xg), ncol = length(yg))

make_contours <- function(xg, yg, z, thresholds, probs, model_name) {
  cl <- contourLines(x = xg, y = yg, z = z, levels = sort(thresholds))
  if (length(cl) == 0L) return(data.frame())
  out <- lapply(seq_along(cl), function(i) {
    k <- which.min(abs(cl[[i]]$level - thresholds))
    data.frame(
      t1 = cl[[i]]$x,
      t2 = cl[[i]]$y,
      model = model_name,
      probability = sprintf("%d%%", round(100 * probs[k])),
      group = paste(model_name, k, i, sep = "_")
    )
  })
  do.call(rbind, out)
}

contours <- rbind(
  make_contours(xg, yg, z_bs2, thr_bs2, hdr_probs, "BS2"),
  make_contours(xg, yg, z_ln, thr_ln, hdr_probs, "Bivariate lognormal")
)
contours$probability <- factor(contours$probability, levels = sprintf("%d%%", round(100 * hdr_probs)))
contours$model <- factor(contours$model, levels = c("BS2", "Bivariate lognormal"))

p <- ggplot(bmd, aes(x = t1, y = t2)) +
  geom_point(size = 2.5, shape = 16) +
  geom_path(data = contours, aes(x = t1, y = t2, group = group), colour = "white", linewidth = 1.35, lineend = "round", inherit.aes = FALSE, show.legend = FALSE) +
  geom_path(data = contours, aes(x = t1, y = t2, group = group, colour = probability, linetype = model, linewidth = model), lineend = "round", inherit.aes = FALSE) +
  geom_text(data = bmd[id == 23], aes(label = id), nudge_x = 0.025, nudge_y = -0.020, size = 3.6) +
  coord_equal(xlim = range(xg), ylim = range(yg), expand = FALSE) +
  labs(x = expression(T[1] ~ "(BMD before)"), y = expression(T[2] ~ "(BMD after)"), colour = "HDR probability", linetype = "Fitted model", linewidth = "Fitted model") +
  scale_linetype_manual(values = c("BS2" = "solid", "Bivariate lognormal" = "longdash")) +
  scale_linewidth_manual(values = c("BS2" = 0.95, "Bivariate lognormal" = 0.75)) +
  guides(colour = guide_legend(order = 2, override.aes = list(linewidth = 1.1)), linetype = guide_legend(order = 1), linewidth = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank())

if ("pdf" %in% formats) ggsave(file.path(fig_dir, "bmd_bs2_vs_lognormal_contours.pdf"), p, width = 7, height = 5.5)
if ("png" %in% formats) ggsave(file.path(fig_dir, "bmd_bs2_vs_lognormal_contours.png"), p, width = 7, height = 5.5, dpi = 600)
cat("Wrote BMD BS2-vs-lognormal contour figure(s) to ", normalizePath(fig_dir), "\n", sep = "")
