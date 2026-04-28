#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# 07_make_figures_eps_publication_v3.R
# Publication-ready EPS figures for the bivariate Birnbaum--Saunders paper.
# Version 3: corrected BMD contour limits and clearer BMD axis labels.
# Run from the root of bs2_parallel_project.
# -----------------------------------------------------------------------------

required <- c("data.table", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "), repos = 'https://cloud.r-project.org')\n",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

message2 <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(pat, "", hit[[1L]])
}

fig_dir <- get_arg("fig-dir", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Output formats. For journal submission the default is EPS only.
# Use --formats=eps,pdf if you also want PDF copies for local preview.
formats <- tolower(trimws(unlist(strsplit(get_arg("formats", "eps"), ",", fixed = TRUE))))
formats <- formats[nzchar(formats)]
if (!all(formats %in% c("eps", "pdf"))) {
  stop("Unsupported --formats value. Use eps, pdf, or eps,pdf.", call. = FALSE)
}

# Number of Monte Carlo samples for the BMD QQ envelope.
B_env <- as.integer(get_arg("B-env", "500"))
seed <- as.integer(get_arg("seed", "202506"))
bmd_file <- get_arg("bmd-file", "data/bmd.csv")


# Axis limits for the BMD contour plot. The default upper limit 1.2 avoids
# clipped contours and leaves the fitted density contours visible beyond
# the observed BMD range.
bmd_contour_min <- as.numeric(get_arg("bmd-contour-min", "0.40"))
bmd_contour_max <- as.numeric(get_arg("bmd-contour-max", "1.20"))
if (!is.finite(bmd_contour_min) || !is.finite(bmd_contour_max) ||
    bmd_contour_min >= bmd_contour_max) {
  stop("Invalid BMD contour limits. Use, for example, --bmd-contour-min=0.40 --bmd-contour-max=1.20",
       call. = FALSE)
}

pick_file <- function(candidates, required = TRUE) {
  candidates <- unique(candidates)
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    if (required) {
      stop("None of these files exists:\n", paste(candidates, collapse = "\n"), call. = FALSE)
    }
    return(NA_character_)
  }
  hit[[1L]]
}

eps_device <- function(filename, width, height, ...) {
  if (capabilities("cairo")) {
    grDevices::cairo_ps(
      filename = filename,
      width = width,
      height = height,
      onefile = FALSE,
      fallback_resolution = 600,
      ...
    )
  } else {
    grDevices::postscript(
      file = filename,
      width = width,
      height = height,
      horizontal = FALSE,
      onefile = FALSE,
      paper = "special",
      ...
    )
  }
}

save_plot <- function(p, filename, width = 7, height = 5, eps = TRUE) {
  # `eps` is kept for backward compatibility. All plots are written in every
  # format listed in --formats, and EPS is the default required format.
  if ("eps" %in% formats) {
    out_eps <- file.path(fig_dir, paste0(filename, ".eps"))
    ggsave(out_eps, plot = p, width = width, height = height, units = "in", device = eps_device)
    message2("Wrote ", out_eps)
  }

  if ("pdf" %in% formats) {
    out_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
    ggsave(out_pdf, plot = p, width = width, height = height, units = "in", device = cairo_pdf)
    message2("Wrote ", out_pdf)
  }
}

clean_stat_labels <- function(x) {
  map <- c(
    LR = "LR",
    Score = "Score",
    Wald = "Wald",
    MM_Wald = "MM-Wald",
    PL_LR = "PL-LR",
    PL_Score = "PL-Score",
    PL_sign = "PL-sign"
  )
  out <- unname(map[x])
  out[is.na(out)] <- gsub("_", "-", x[is.na(out)])
  out
}

hypothesis_label <- function(x) {
  map <- c(H01 = "H01", H02 = "H02", H03 = "H03", H04 = "H04")
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

rho_label <- function(x) paste0("rho = ", format(x, trim = TRUE, scientific = FALSE))
n_label <- function(x) paste0("n = ", x)

# Ordered factor helpers for publication facets and discrete sample-size axes.
# These avoid alphabetical facet ordering such as n = 120, n = 30, n = 60.
n_factor <- function(x, levels = c(24, 30, 60, 120)) {
  factor(x, levels = levels, labels = paste0("n = ", levels))
}

n_tick_factor <- function(x, levels = c(24, 30, 60, 120)) {
  factor(x, levels = levels, labels = as.character(levels))
}

rho_factor <- function(x, levels = c(0, 0.6)) {
  factor(x, levels = levels, labels = paste0("rho = ", format(levels, trim = TRUE, scientific = FALSE)))
}

alt_factor <- function(x) {
  factor(
    x,
    levels = c("shape_only", "scale_only", "shape_and_scale"),
    labels = c("Shape only", "Scale only", "Shape and scale")
  )
}

alt_label <- function(x) {
  map <- c(
    scale_only = "Scale only",
    shape_only = "Shape only",
    shape_and_scale = "Shape and scale"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

dgp_label_clean <- function(x) {
  map <- c(
    bs2 = "BS2",
    contaminated_bs2 = "Contaminated BS2",
    lognormal = "Lognormal",
    BS2 = "BS2",
    `Contaminated BS2` = "Contaminated BS2",
    Lognormal = "Lognormal"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

# Black-and-white line/shape scheme, suitable for EPS and printed articles.
stat_levels_all <- c("LR", "Score", "Wald", "MM-Wald", "PL-LR", "PL-Score", "PL-sign")
stat_linetypes <- c(
  "LR" = "solid",
  "Score" = "dotdash",
  "Wald" = "longdash",
  "MM-Wald" = "dashed",
  "PL-LR" = "twodash",
  "PL-Score" = "solid",
  "PL-sign" = "dotted"
)
stat_shapes <- c(
  "LR" = 16,
  "Score" = 8,
  "Wald" = 3,
  "MM-Wald" = 17,
  "PL-LR" = 15,
  "PL-Score" = 4,
  "PL-sign" = 7
)

base_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "grey55"),
      strip.text = element_text(size = base_size - 1),
      plot.title = element_blank(),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1)
    )
}

apply_stat_scales <- function(levels = stat_levels_all) {
  list(
    scale_linetype_manual(values = stat_linetypes[levels], breaks = levels, drop = FALSE),
    scale_shape_manual(values = stat_shapes[levels], breaks = levels, drop = FALSE)
  )
}

# -----------------------------------------------------------------------------
# BS2 utilities for BMD figures
# -----------------------------------------------------------------------------

source_core_if_available <- function() {
  if (file.exists("R/bs2_core.R")) {
    source("R/bs2_core.R")
    return(invisible(TRUE))
  }
  warning("R/bs2_core.R not found. BMD fit functions must already be loaded.")
  invisible(FALSE)
}

source_core_if_available()

# Fallback definitions in case the core file uses different helper names.
if (!exists("a_bs", mode = "function")) {
  a_bs <- function(t, alpha, beta) {
    (sqrt(t / beta) - sqrt(beta / t)) / alpha
  }
}

if (!exists("A_bs", mode = "function")) {
  A_bs <- function(t, alpha, beta) {
    (t + beta) / (2 * alpha * sqrt(beta) * t^(3 / 2))
  }
}

if (!exists("ell_bs2", mode = "function")) {
  ell_bs2 <- function(t1, t2, alpha, beta, rho) {
    if (any(alpha <= 0) || any(beta <= 0) || abs(rho) >= 1) return(-Inf)
    z1 <- a_bs(t1, alpha[1], beta[1])
    z2 <- a_bs(t2, alpha[2], beta[2])
    q <- (z1^2 + z2^2 - 2 * rho * z1 * z2) / (1 - rho^2)
    sum(-log(2 * pi) - 0.5 * log(1 - rho^2) - 0.5 * q +
          log(A_bs(t1, alpha[1], beta[1])) + log(A_bs(t2, alpha[2], beta[2])))
  }
}

if (!exists("dbs2", mode = "function")) {
  dbs2 <- function(t1, t2, alpha, beta, rho, log = FALSE) {
    z1 <- a_bs(t1, alpha[1], beta[1])
    z2 <- a_bs(t2, alpha[2], beta[2])
    q <- (z1^2 + z2^2 - 2 * rho * z1 * z2) / (1 - rho^2)
    val <- -log(2 * pi) - 0.5 * log(1 - rho^2) - 0.5 * q +
      log(A_bs(t1, alpha[1], beta[1])) + log(A_bs(t2, alpha[2], beta[2]))
    if (log) val else exp(val)
  }
}

bs2_mahalanobis <- function(t1, t2, fit) {
  z1 <- a_bs(t1, fit$alpha[1], fit$beta[1])
  z2 <- a_bs(t2, fit$alpha[2], fit$beta[2])
  rho <- fit$rho
  (z1^2 + z2^2 - 2 * rho * z1 * z2) / (1 - rho^2)
}

bmd_data <- function(file = bmd_file) {
  if (!file.exists(file)) {
    stop("BMD data file not found: ", file,
         "\nRun from the repository root or pass --bmd-file=path/to/bmd.csv.",
         call. = FALSE)
  }
  dat <- data.table::fread(file)
  required <- c("id", "t1", "t2")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("BMD data file is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dat[, .(id = as.integer(id), t1 = as.numeric(t1), t2 = as.numeric(t2))]
}

fit_bmd_bs2 <- function(dat) {
  if (!exists("fit_bs2_mm", mode = "function") || !exists("fit_bs2_profile", mode = "function")) {
    stop("fit_bs2_mm() and fit_bs2_profile() are required. Source R/bs2_core.R first.", call. = FALSE)
  }
  mm <- fit_bs2_mm(dat$t1, dat$t2)
  fit <- fit_bs2_profile(dat$t1, dat$t2, beta_start = mm$beta)
  fit
}

make_bmd_envelope <- function(dat, fit, B = 500, seed = 202506) {
  if (!exists("rbs2", mode = "function")) stop("rbs2() is required. Source R/bs2_core.R first.", call. = FALSE)
  set.seed(seed)
  n <- nrow(dat)
  theo <- qchisq(ppoints(n), df = 2)
  obs <- sort(bs2_mahalanobis(dat$t1, dat$t2, fit))
  env <- matrix(NA_real_, nrow = B, ncol = n)

  for (b in seq_len(B)) {
    xb <- rbs2(n, alpha = fit$alpha, beta = fit$beta, rho = fit$rho)
    xb <- as.data.table(xb)
    setnames(xb, c("t1", "t2"))
    fb <- tryCatch({
      mmb <- fit_bs2_mm(xb$t1, xb$t2)
      fit_bs2_profile(xb$t1, xb$t2, beta_start = mmb$beta)
    }, error = function(e) NULL)
    if (!is.null(fb)) env[b, ] <- sort(bs2_mahalanobis(xb$t1, xb$t2, fb))
  }

  qfun <- function(prob) apply(env, 2, quantile, probs = prob, na.rm = TRUE, type = 8)
  env_dt <- data.table(
    theo = theo,
    lower = qfun(0.05),
    median = qfun(0.50),
    upper = qfun(0.95),
    obs = obs
  )
  env_lines <- melt(
    env_dt,
    id.vars = "theo",
    measure.vars = c("lower", "median", "upper"),
    variable.name = "envelope",
    value.name = "distance"
  )
  env_lines[, envelope := factor(
    envelope,
    levels = c("lower", "median", "upper"),
    labels = c("5% envelope", "Median envelope", "95% envelope")
  )]

  ggplot() +
    geom_line(data = env_lines, aes(x = theo, y = distance, linetype = envelope),
              linewidth = 0.45) +
    geom_point(data = env_dt, aes(x = theo, y = obs, shape = "Observed"), size = 1.8) +
    scale_linetype_manual(
      values = c("5% envelope" = "dashed", "Median envelope" = "solid", "95% envelope" = "dashed"),
      name = NULL
    ) +
    scale_shape_manual(values = c("Observed" = 16), name = NULL) +
    labs(
      x = "Chi-square(2) quantiles",
      y = "Ordered fitted Mahalanobis distances"
    ) +
    guides(
      linetype = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    base_theme(12) +
    theme(legend.position = "bottom")
}

make_bmd_contours <- function(dat, fit, lower = bmd_contour_min, upper = bmd_contour_max) {
  # The grid is intentionally wider than the observed BMD range. In the
  # manuscript figure this prevents the fitted contours from being truncated
  # at the upper/right plotting borders.
  t1_grid <- seq(lower, upper, length.out = 260)
  t2_grid <- seq(lower, upper, length.out = 260)
  gr <- CJ(t1 = t1_grid, t2 = t2_grid)
  gr[, dens := dbs2(t1, t2, alpha = fit$alpha, beta = fit$beta, rho = fit$rho)]

  ggplot() +
    geom_contour(data = gr, aes(x = t1, y = t2, z = dens),
                 bins = 8, linewidth = 0.45, colour = "grey25") +
    geom_point(data = dat, aes(x = t1, y = t2), size = 2) +
    geom_text(data = dat[id == 23], aes(x = t1, y = t2, label = "23"),
              nudge_x = 0.025, nudge_y = 0.025, size = 3.5) +
    labs(
      x = "BMD before study (T1)",
      y = "BMD after study (T2)"
    ) +
    coord_equal(xlim = c(lower, upper), ylim = c(lower, upper), expand = FALSE) +
    scale_x_continuous(breaks = seq(lower, upper, by = 0.2)) +
    scale_y_continuous(breaks = seq(lower, upper, by = 0.2)) +
    base_theme(12) +
    theme(legend.position = "none")
}


make_bmd_case_deletion <- function(dat, fit) {
  theta0 <- c(fit$alpha, fit$beta, fit$rho)
  names(theta0) <- c("alpha1", "alpha2", "beta1", "beta2", "rho")
  out <- vector("list", nrow(dat))

  for (i in seq_len(nrow(dat))) {
    dd <- dat[-i]
    fd <- tryCatch(fit_bmd_bs2(dd), error = function(e) NULL)
    if (is.null(fd)) {
      out[[i]] <- data.table(id = i, ok = FALSE, rel_norm = NA_real_)
    } else {
      th <- c(fd$alpha, fd$beta, fd$rho)
      delta <- th - theta0
      rel_norm <- sqrt(sum((delta / pmax(abs(theta0), 1e-8))^2))
      out[[i]] <- data.table(
        id = i,
        ok = TRUE,
        delta_alpha1 = delta[1],
        delta_alpha2 = delta[2],
        delta_beta1 = delta[3],
        delta_beta2 = delta[4],
        delta_rho = delta[5],
        rel_norm = rel_norm
      )
    }
  }

  infl <- rbindlist(out, fill = TRUE)
  fwrite(infl, file.path(fig_dir, "bmd_case_deletion_values.csv"))

  y_upper <- max(0.52, max(infl$rel_norm, na.rm = TRUE) + 0.06)

  ggplot(infl, aes(x = id, y = rel_norm)) +
    geom_segment(aes(xend = id, y = 0, yend = rel_norm), linewidth = 0.35) +
    geom_point(aes(shape = id == 23), size = 2.4) +
    geom_text(data = infl[id == 23], aes(label = "23"), nudge_y = 0.035, size = 3.5) +
    scale_x_continuous(breaks = seq(1, nrow(dat), by = 2)) +
    scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0, 0.02))) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    labs(
      x = "Deleted observation",
      y = "Relative perturbation norm"
    ) +
    base_theme(12) +
    theme(legend.position = "none")
}

# -----------------------------------------------------------------------------
# Figure 1: BMD envelope, BMD contours, BMD influence
# -----------------------------------------------------------------------------

message2("Generating BMD diagnostic figures...")
bmd <- bmd_data(bmd_file)
fit_bmd <- fit_bmd_bs2(bmd)
message2("BMD fit: alpha=", paste(round(fit_bmd$alpha, 4), collapse = ","),
         "; beta=", paste(round(fit_bmd$beta, 4), collapse = ","),
         "; rho=", round(fit_bmd$rho, 4))

p_env <- make_bmd_envelope(bmd, fit_bmd, B = B_env, seed = seed)
save_plot(p_env, "envelope", width = 4.8, height = 4.0, eps = TRUE)

p_cont <- make_bmd_contours(bmd, fit_bmd)
save_plot(p_cont, "contorno", width = 4.8, height = 4.0, eps = TRUE)

p_infl <- make_bmd_case_deletion(bmd, fit_bmd)
save_plot(p_infl, "bmd_case_deletion", width = 6.8, height = 4.2, eps = FALSE)

# -----------------------------------------------------------------------------
# Figure 2: empirical size
# -----------------------------------------------------------------------------

message2("Generating empirical size figure...")
size_file <- pick_file(c(
  "results/tables_score_fixed/test_rejection_summary.csv",
  "results/tables_score_fixed/test_rejection_summary_score_fixed_clean.csv",
  "results/tables/clean/test_rejection_summary_clean.csv",
  "test_rejection_summary_score_fixed_clean.csv",
  "test_rejection_summary_clean.csv"
), required = FALSE)

if (!is.na(size_file)) {
  size_dt <- fread(size_file)
  size_dt <- size_dt[!(hypothesis != "H04" & grepl("^PL", statistic))]
  size_dt[, statistic := factor(statistic,
                                levels = c("LR", "Score", "Wald", "MM_Wald", "PL_LR", "PL_Score", "PL_sign"))]

  # Average over nuisance rho configurations for H01--H03, while H04 has one null curve.
  size_plot_dt <- size_dt[!is.na(statistic), .(
    rejection_rate = mean(rejection_rate, na.rm = TRUE),
    R_eff = sum(R_eff, na.rm = TRUE)
  ), by = .(hypothesis, n, statistic)]
  size_plot_dt[, stat_label := factor(clean_stat_labels(as.character(statistic)), levels = stat_levels_all)]
  size_plot_dt[, hypothesis_label := factor(hypothesis_label(hypothesis), levels = c("H01", "H02", "H03", "H04"))]
  # Use a categorical x-axis so that the close values 24 and 30 do not overlap.
  size_plot_dt[, n_fac := n_tick_factor(n, levels = c(24, 30, 60, 120))]

  band_low <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / 3000)
  band_high <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / 3000)

  size_levels <- stat_levels_all[stat_levels_all %in% levels(droplevels(size_plot_dt$stat_label))]
  p_size <- ggplot(size_plot_dt, aes(x = n_fac, y = rejection_rate,
                                     linetype = stat_label, shape = stat_label,
                                     group = stat_label)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = band_low, ymax = band_high,
             fill = "grey90", colour = NA) +
    geom_hline(yintercept = 0.05, linewidth = 0.35) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.8) +
    facet_wrap(~ hypothesis_label, nrow = 2) +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Sample size",
      y = "Empirical rejection rate",
      linetype = "Statistic",
      shape = "Statistic"
    ) +
    apply_stat_scales(size_levels) +
    base_theme(12)

  save_plot(p_size, "empirical_size", width = 7.4, height = 5.5)
} else {
  warning("Size results file not found; empirical_size.eps was not generated.")
}

# -----------------------------------------------------------------------------
# Figure 3: empirical power
# -----------------------------------------------------------------------------

message2("Generating empirical power figures...")
power_file <- pick_file(c(
  "results/tables_power_score_fixed/power_rejection_summary.csv",
  "results/tables_power_score_fixed/power_rejection_summary_reviewed.csv",
  "power_rejection_summary_reviewed.csv"
), required = FALSE)

if (!is.na(power_file)) {
  pow <- fread(power_file)
  pow[, stat_label := factor(clean_stat_labels(statistic), levels = stat_levels_all)]
  # Ordered factors control panel order in all power plots.
  pow[, rho_fac := rho_factor(rho_true, levels = c(0, 0.6))]
  pow[, n_fac := n_factor(n, levels = c(24, 30, 60, 120))]
  if ("alt_id" %in% names(pow)) pow[, alt_fac := alt_factor(alt_id)]

  classical <- c("LR", "Score", "Wald", "MM_Wald")
  classical_labels <- clean_stat_labels(classical)
  h04_stats <- c("LR", "Score", "Wald", "MM_Wald", "PL_LR", "PL_Score", "PL_sign")
  h04_labels <- clean_stat_labels(h04_stats)

  # H01: power against alpha ratio alternatives.
  p_h01 <- ggplot(pow[hypothesis == "H01" & true_null == FALSE & statistic %in% classical],
                  aes(x = alpha_ratio, y = rejection_rate,
                      linetype = stat_label, shape = stat_label,
                      group = stat_label)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.7) +
    facet_grid(rho_fac ~ n_fac, drop = TRUE) +
    scale_x_continuous(breaks = c(1.25, 1.5, 2.0)) +
    labs(x = expression(alpha[2]/alpha[1]), y = "Empirical power",
         linetype = "Statistic", shape = "Statistic") +
    apply_stat_scales(classical_labels) +
    base_theme(11)
  save_plot(p_h01, "power_H01", width = 7.6, height = 4.8)

  # H02: power against beta ratio alternatives.
  p_h02 <- ggplot(pow[hypothesis == "H02" & true_null == FALSE & statistic %in% classical],
                  aes(x = beta_ratio, y = rejection_rate,
                      linetype = stat_label, shape = stat_label,
                      group = stat_label)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.7) +
    facet_grid(rho_fac ~ n_fac, drop = TRUE) +
    scale_x_continuous(breaks = c(1.25, 1.5, 2.0)) +
    labs(x = expression(beta[2]/beta[1]), y = "Empirical power",
         linetype = "Statistic", shape = "Statistic") +
    apply_stat_scales(classical_labels) +
    base_theme(11)
  save_plot(p_h02, "power_H02", width = 7.6, height = 4.8)

  # H03: power for joint alternatives.
  p_h03 <- ggplot(pow[hypothesis == "H03" & true_null == FALSE & statistic %in% classical],
                  aes(x = n_tick_factor(n, levels = c(30, 60, 120)), y = rejection_rate,
                      linetype = stat_label, shape = stat_label,
                      group = stat_label)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.7) +
    facet_grid(alt_fac ~ rho_fac, drop = TRUE) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Sample size", y = "Empirical power",
         linetype = "Statistic", shape = "Statistic") +
    apply_stat_scales(classical_labels) +
    base_theme(11)
  save_plot(p_h03, "power_H03", width = 7.8, height = 6.4)

  # H04: power for independence, including PL procedures.
  p_h04 <- ggplot(pow[hypothesis == "H04" & true_null == FALSE & statistic %in% h04_stats],
                  aes(x = rho_true, y = rejection_rate,
                      linetype = stat_label, shape = stat_label,
                      group = stat_label)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.7) +
    facet_wrap(~ n_fac, nrow = 1, drop = TRUE) +
    scale_x_continuous(breaks = c(-0.6, -0.3, 0.3, 0.6)) +
    labs(x = expression(rho), y = "Empirical power",
         linetype = "Statistic", shape = "Statistic") +
    apply_stat_scales(h04_labels) +
    base_theme(11)
  save_plot(p_h04, "power_H04", width = 8.8, height = 4.2)
} else {
  warning("Power results file not found; EPS power figures were not generated.")
}

# -----------------------------------------------------------------------------
# Figure 4: GOF/influence simulation
# -----------------------------------------------------------------------------

message2("Generating GOF/influence figure...")
gof_file <- pick_file(c(
  "results/tables_gof_influence/gof_influence_manuscript_summary.csv",
  "results/tables_gof_influence/gof_influence_summary_clean.csv",
  "gof_influence_summary_clean.csv"
), required = FALSE)

if (!is.na(gof_file)) {
  gof <- fread(gof_file)
  if (!"dgp_label" %in% names(gof)) gof[, dgp_label := dgp]
  gof[, dgp_label := factor(dgp_label_clean(dgp_label),
                            levels = c("BS2", "Contaminated BS2", "Lognormal"))]
  gof_long <- melt(
    gof,
    id.vars = c("dgp", "dgp_label", "n"),
    measure.vars = c("rejection_rate", "median_max_rel_norm"),
    variable.name = "quantity",
    value.name = "value"
  )
  gof_long[, quantity := fifelse(quantity == "rejection_rate",
                                 "GOF rejection rate",
                                 "Median max. relative perturbation")]

  gof_hline <- data.table(quantity = "GOF rejection rate", yintercept = 0.05)

  p_gof <- ggplot(gof_long, aes(x = n, y = value,
                                linetype = dgp_label, shape = dgp_label)) +
    geom_hline(data = gof_hline, aes(yintercept = yintercept),
               inherit.aes = FALSE, linewidth = 0.35) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2.0) +
    facet_wrap(~ quantity, scales = "free_y", nrow = 1) +
    scale_x_continuous(breaks = c(24, 50, 100)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Sample size", y = NULL,
         linetype = "DGP", shape = "DGP") +
    scale_linetype_manual(values = c("BS2" = "solid", "Contaminated BS2" = "dashed", "Lognormal" = "dotdash"), drop = FALSE) +
    scale_shape_manual(values = c("BS2" = 16, "Contaminated BS2" = 17, "Lognormal" = 15), drop = FALSE) +
    base_theme(12)

  save_plot(p_gof, "gof_influence_sim", width = 8.4, height = 4.2)
} else {
  warning("GOF/influence results file not found; gof_influence_sim.eps was not generated.")
}

message2("Done. EPS figures written to: ", normalizePath(fig_dir))
