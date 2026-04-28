# Core helper functions for the bivariate Birnbaum-Saunders (BS2) simulation code.
# These functions use only base R and stats.
# Source this file before sourcing fit_bs2_restricted.R, compute_bs2_tests.R,
# boot_lr_h0.R, or fit_bs2_pl_crossfit.R.

.bs2_check_data <- function(t1, t2) {
  if (length(t1) != length(t2)) stop("t1 and t2 must have the same length.")
  if (length(t1) < 2L) stop("At least two observations are required.")
  if (any(!is.finite(t1)) || any(!is.finite(t2))) stop("Data must be finite.")
  if (any(t1 <= 0) || any(t2 <= 0)) stop("BS2 data must be strictly positive.")
  invisible(TRUE)
}

.clip_rho <- function(rho, eps = 1e-8) {
  rho <- as.numeric(rho)[1L]
  if (!is.finite(rho)) stop("rho must be finite.")
  max(min(rho, 1 - eps), -1 + eps)
}

a_bs <- function(t, alpha, beta) {
  (sqrt(t / beta) - sqrt(beta / t)) / alpha
}

a_bs_unit_alpha <- function(t, beta) {
  sqrt(t / beta) - sqrt(beta / t)
}

A_bs <- function(t, alpha, beta) {
  (t + beta) / (2 * alpha * sqrt(beta) * t^(3 / 2))
}

rho_from_lambda <- function(lambda) {
  lambda <- as.numeric(lambda)
  out <- lambda / (1 + sqrt(1 + lambda^2))
  out[is.infinite(lambda) & lambda > 0] <- 1
  out[is.infinite(lambda) & lambda < 0] <- -1
  out
}

lambda_from_rho <- function(rho) {
  rho <- .clip_rho(rho, eps = 1e-12)
  2 * rho / (1 - rho^2)
}

log1pexp <- function(z) {
  ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
}

.fit_ok <- function(x) {
  !is.null(x) && !inherits(x, "try-error")
}

.safe_eval <- function(expr) {
  tryCatch(expr, error = function(e) structure(list(error = conditionMessage(e)), class = "try-error"))
}

.safe_solve <- function(M, tol = sqrt(.Machine$double.eps)) {
  M <- as.matrix(M)
  if (!all(is.finite(M))) stop("Matrix contains non-finite values.")
  M <- (M + t(M)) / 2
  s <- svd(M)
  d <- s$d
  if (length(d) == 0L || max(abs(d)) == 0) stop("Matrix is numerically singular.")
  keep <- d > tol * max(d)
  if (!any(keep)) stop("Matrix is numerically singular.")
  s$v[, keep, drop = FALSE] %*% (diag(1 / d[keep], nrow = sum(keep)) %*% t(s$u[, keep, drop = FALSE]))
}

.fit_to_theta <- function(fit) {
  c(as.numeric(fit$alpha), as.numeric(fit$beta), as.numeric(fit$rho))
}

.restriction_matrix <- function(hypothesis, rho0 = 0) {
  h <- toupper(gsub("[^A-Z0-9]", "", hypothesis))
  if (h %in% c("H1", "H01")) h <- "H01"
  if (h %in% c("H2", "H02")) h <- "H02"
  if (h %in% c("H3", "H03")) h <- "H03"
  if (h %in% c("H4", "H04")) h <- "H04"

  if (h == "H01") {
    A <- matrix(c(1, -1, 0, 0, 0), nrow = 1)
    q <- 0
  } else if (h == "H02") {
    A <- matrix(c(0, 0, 1, -1, 0), nrow = 1)
    q <- 0
  } else if (h == "H03") {
    A <- rbind(c(1, -1, 0, 0, 0), c(0, 0, 1, -1, 0))
    q <- c(0, 0)
  } else if (h == "H04") {
    A <- matrix(c(0, 0, 0, 0, 1), nrow = 1)
    q <- rho0
  } else {
    stop("Unknown hypothesis. Use H01, H02, H03, or H04.")
  }
  list(A = A, q = as.numeric(q), df = nrow(A), hypothesis = h)
}

rbs2 <- function(n, alpha, beta, rho) {
  if (length(alpha) != 2L || length(beta) != 2L) stop("alpha and beta must have length 2.")
  if (any(alpha <= 0) || any(beta <= 0)) stop("alpha and beta must be positive.")
  rho <- .clip_rho(rho, eps = 1e-12)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  C <- chol(Sigma)
  z <- matrix(rnorm(2 * n), nrow = n, ncol = 2) %*% C

  t1 <- beta[1] * (
    alpha[1] * z[, 1] / 2 + sqrt((alpha[1] * z[, 1] / 2)^2 + 1)
  )^2
  t2 <- beta[2] * (
    alpha[2] * z[, 2] / 2 + sqrt((alpha[2] * z[, 2] / 2)^2 + 1)
  )^2
  cbind(t1, t2)
}

ell_bs2 <- function(t1, t2, alpha, beta, rho) {
  if (length(alpha) != 2L || length(beta) != 2L) return(-Inf)
  if (any(!is.finite(alpha)) || any(!is.finite(beta)) || !is.finite(rho)) return(-Inf)
  if (any(alpha <= 0) || any(beta <= 0) || abs(rho) >= 1) return(-Inf)
  if (any(t1 <= 0) || any(t2 <= 0)) return(-Inf)

  z1 <- a_bs(t1, alpha[1], beta[1])
  z2 <- a_bs(t2, alpha[2], beta[2])
  A1 <- A_bs(t1, alpha[1], beta[1])
  A2 <- A_bs(t2, alpha[2], beta[2])
  if (any(A1 <= 0) || any(A2 <= 0)) return(-Inf)

  den <- 1 - rho^2
  q <- (z1^2 + z2^2 - 2 * rho * z1 * z2) / den

  val <- sum(-log(2 * pi) - 0.5 * log(den) - 0.5 * q + log(A1) + log(A2))
  if (!is.finite(val)) -Inf else val
}

fit_bs2_mm <- function(t1, t2, rho_eps = 1e-8) {
  .bs2_check_data(t1, t2)

  S1 <- mean(t1)
  S2 <- mean(t2)
  R1 <- 1 / mean(1 / t1)
  R2 <- 1 / mean(1 / t2)

  alpha1 <- sqrt(max(0, 2 * (sqrt(S1 / R1) - 1)))
  alpha2 <- sqrt(max(0, 2 * (sqrt(S2 / R2) - 1)))
  beta1 <- sqrt(S1 * R1)
  beta2 <- sqrt(S2 * R2)

  q1 <- a_bs_unit_alpha(t1, beta1)
  q2 <- a_bs_unit_alpha(t2, beta2)
  den <- sqrt(sum(q1^2) * sum(q2^2))
  rho <- if (den > 0) sum(q1 * q2) / den else 0
  rho <- .clip_rho(rho, eps = rho_eps)

  list(alpha = c(alpha1, alpha2), beta = c(beta1, beta2), rho = rho)
}

.bs2_profile_unrestricted_from_beta <- function(t1, t2, beta, rho_eps = 1e-8) {
  if (length(beta) != 2L || any(!is.finite(beta)) || any(beta <= 0)) return(NULL)
  q1 <- a_bs_unit_alpha(t1, beta[1])
  q2 <- a_bs_unit_alpha(t2, beta[2])
  s11 <- sum(q1^2)
  s22 <- sum(q2^2)
  if (s11 <= 0 || s22 <= 0) return(NULL)
  alpha <- c(sqrt(s11 / length(t1)), sqrt(s22 / length(t1)))
  rho <- sum(q1 * q2) / sqrt(s11 * s22)
  rho <- .clip_rho(rho, eps = rho_eps)
  logLik <- ell_bs2(t1, t2, alpha = alpha, beta = beta, rho = rho)
  list(alpha = alpha, beta = beta, rho = rho, logLik = logLik)
}

.optim_best <- function(par, fn, methods = c("BFGS", "Nelder-Mead"), control = list(maxit = 1000)) {
  best <- NULL
  for (method in methods) {
    opt <- try(optim(par = par, fn = fn, method = method, control = control), silent = TRUE)
    if (inherits(opt, "try-error")) next
    if (is.null(best) || is.finite(opt$value) && opt$value < best$value) best <- opt
  }
  if (is.null(best)) stop("All optim() attempts failed.")
  best
}

fit_bs2_profile <- function(t1, t2, beta_start = NULL, rho_eps = 1e-8,
                            optim_methods = c("BFGS", "Nelder-Mead"),
                            control = list(maxit = 1000)) {
  .bs2_check_data(t1, t2)
  if (is.null(beta_start)) beta_start <- fit_bs2_mm(t1, t2)$beta
  if (length(beta_start) != 2L) stop("beta_start must have length 2.")
  beta_start <- pmax(as.numeric(beta_start), .Machine$double.eps)

  loglik_prof <- function(eta_beta) {
    beta <- exp(eta_beta)
    fit <- .bs2_profile_unrestricted_from_beta(t1, t2, beta = beta, rho_eps = rho_eps)
    if (is.null(fit)) return(-Inf)
    fit$logLik
  }

  objective <- function(eta_beta) {
    val <- loglik_prof(eta_beta)
    if (!is.finite(val)) return(.Machine$double.xmax / 100)
    -val
  }

  opt <- .optim_best(log(beta_start), objective, methods = optim_methods, control = control)
  beta_hat <- exp(opt$par)
  fit <- .bs2_profile_unrestricted_from_beta(t1, t2, beta = beta_hat, rho_eps = rho_eps)

  fit$convergence <- opt$convergence
  fit$message <- opt$message
  fit$objective <- opt$value
  fit$theta <- .fit_to_theta(fit)
  fit
}

# One-dimensional conditional logistic slope fit used by fit_bs2_pl() and cross-fitting.
.fit_conditional_logistic_slope <- function(x, c, penalized = FALSE, lambda_max = 50) {
  if (length(x) != length(c)) stop("x and c must have the same length.")
  if (any(!is.finite(x)) || any(x < 0)) stop("x must be finite and non-negative.")
  if (any(!(c %in% c(0, 1)))) stop("c must be a 0/1 vector.")

  ell <- function(lambda) {
    eta <- lambda * x
    sum(c * eta - log1pexp(eta))
  }
  info <- function(lambda) {
    pi_hat <- plogis(lambda * x)
    sum(x^2 * pi_hat * (1 - pi_hat))
  }
  crit <- function(lambda) {
    val <- ell(lambda)
    if (penalized) {
      j <- info(lambda)
      if (!is.finite(j) || j <= 0) return(-Inf)
      val <- val + 0.5 * log(j)
    }
    val
  }

  separation <- all(c == 1) || all(c == 0)
  if (separation && !penalized) {
    lambda_hat <- if (all(c == 1)) Inf else -Inf
    rho_hat <- if (lambda_hat > 0) 1 else -1
    logLik_hat <- ell(sign(lambda_hat) * lambda_max)
    return(list(lambda = lambda_hat, rho = rho_hat, logLik = logLik_hat,
                separation = TRUE, boundary = TRUE, convergence = NA_integer_))
  }

  opt <- optimize(function(lambda) -crit(lambda), interval = c(-lambda_max, lambda_max))
  lambda_hat <- opt$minimum
  rho_hat <- rho_from_lambda(lambda_hat)
  boundary <- abs(lambda_hat) > 0.98 * lambda_max

  list(lambda = lambda_hat, rho = rho_hat, logLik = ell(lambda_hat),
       separation = separation, boundary = boundary, convergence = 0L)
}

fit_bs2_pl <- function(t1, t2, alpha_hat, beta_hat, penalized = FALSE, lambda_max = 50) {
  .bs2_check_data(t1, t2)
  y1 <- a_bs(t1, alpha_hat[1], beta_hat[1])
  y2 <- a_bs(t2, alpha_hat[2], beta_hat[2])
  x <- abs(y1 * y2)
  cind <- as.numeric(y1 * y2 > 0)

  fit <- .fit_conditional_logistic_slope(x, cind, penalized = penalized, lambda_max = lambda_max)
  lambda_hat <- fit$lambda
  rho_hat <- fit$rho

  ell0 <- sum(cind * 0 - log1pexp(0 * x))
  lr0 <- if (is.finite(lambda_hat)) 2 * (fit$logLik - ell0) else NA_real_

  pi_hat <- if (is.finite(lambda_hat)) plogis(lambda_hat * x) else rep(ifelse(lambda_hat > 0, 1, 0), length(x))
  A_hat <- mean(x^2 * pi_hat * (1 - pi_hat))
  B_hat <- mean(x^2 * (cind - pi_hat)^2)

  se_lambda <- if (is.finite(A_hat) && A_hat > 0) sqrt((1 / length(x)) * B_hat / A_hat^2) else NA_real_
  gprime <- if (abs(rho_hat) < 1) (1 - rho_hat^2)^2 / (2 * (1 + rho_hat^2)) else NA_real_
  se_rho <- if (is.finite(gprime) && is.finite(se_lambda)) abs(gprime) * se_lambda else NA_real_

  score0 <- if (sum(x^2) > 0) (sum(x * (2 * cind - 1))^2) / sum(x^2) else NA_real_

  list(lambda = lambda_hat, rho = rho_hat, se_lambda = se_lambda, se_rho = se_rho,
       logLik = fit$logLik, lr0 = lr0, score0 = score0, x = x, c = cind,
       separation = fit$separation, boundary = fit$boundary, convergence = fit$convergence)
}

calibrate_pl_score_fast <- function(x, c, B = 9999L, chunk = 2000L) {
  s_obs <- (sum(x * (2 * c - 1))^2) / sum(x^2)
  ge <- 0L
  done <- 0L
  n <- length(x)
  denom <- sum(x^2)
  while (done < B) {
    b <- min(chunk, B - done)
    signs <- matrix(sample(c(-1, 1), size = b * n, replace = TRUE), nrow = b, ncol = n)
    s_star <- as.vector(signs %*% x)^2 / denom
    ge <- ge + sum(s_star >= s_obs)
    done <- done + b
  }
  (1 + ge) / (B + 1)
}

.bs2_theta_valid <- function(theta) {
  length(theta) == 5L && all(is.finite(theta)) && all(theta[1:4] > 0) && abs(theta[5]) < 1
}

.bs2_loglik_theta <- function(theta, t1, t2) {
  if (!.bs2_theta_valid(theta)) return(-Inf)
  ell_bs2(t1, t2, alpha = theta[1:2], beta = theta[3:4], rho = theta[5])
}

.bs2_fd_steps <- function(theta, rel_step = 1e-5) {
  h <- rel_step * pmax(abs(theta), 1)
  h[1:4] <- pmin(h[1:4], theta[1:4] / 4)
  rho_room <- max(1e-8, (1 - abs(theta[5])) / 4)
  h[5] <- min(h[5], rho_room)
  pmax(h, 1e-7)
}

.num_grad <- function(f, x, h = NULL) {
  x <- as.numeric(x)
  p <- length(x)
  if (is.null(h)) h <- .bs2_fd_steps(x)
  g <- numeric(p)
  for (j in seq_len(p)) {
    xp <- x; xm <- x
    xp[j] <- xp[j] + h[j]
    xm[j] <- xm[j] - h[j]
    fp <- f(xp); fm <- f(xm)
    if (!is.finite(fp) || !is.finite(fm)) {
      h2 <- h[j] / 10
      xp <- x; xm <- x
      xp[j] <- xp[j] + h2
      xm[j] <- xm[j] - h2
      fp <- f(xp); fm <- f(xm)
      h[j] <- h2
    }
    g[j] <- (fp - fm) / (2 * h[j])
  }
  g
}

.num_hessian <- function(f, x, h = NULL) {
  x <- as.numeric(x)
  p <- length(x)
  if (is.null(h)) h <- .bs2_fd_steps(x)
  H <- matrix(NA_real_, p, p)
  f0 <- f(x)

  for (i in seq_len(p)) {
    xp <- x; xm <- x
    xp[i] <- xp[i] + h[i]
    xm[i] <- xm[i] - h[i]
    fp <- f(xp); fm <- f(xm)
    H[i, i] <- (fp - 2 * f0 + fm) / (h[i]^2)
  }

  if (p >= 2L) {
    for (i in seq_len(p - 1L)) {
      for (j in (i + 1L):p) {
        xpp <- x; xpm <- x; xmp <- x; xmm <- x
        xpp[i] <- xpp[i] + h[i]; xpp[j] <- xpp[j] + h[j]
        xpm[i] <- xpm[i] + h[i]; xpm[j] <- xpm[j] - h[j]
        xmp[i] <- xmp[i] - h[i]; xmp[j] <- xmp[j] + h[j]
        xmm[i] <- xmm[i] - h[i]; xmm[j] <- xmm[j] - h[j]
        H[i, j] <- (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h[i] * h[j])
        H[j, i] <- H[i, j]
      }
    }
  }
  (H + t(H)) / 2
}

score_bs2 <- function(t1, t2, theta, rel_step = 1e-5) {
  .bs2_check_data(t1, t2)
  if (!.bs2_theta_valid(theta)) stop("Invalid theta.")
  f <- function(th) .bs2_loglik_theta(th, t1, t2)
  .num_grad(f, theta, h = .bs2_fd_steps(theta, rel_step = rel_step))
}

obs_info_bs2 <- function(t1, t2, theta, rel_step = 1e-5, ridge = 1e-8) {
  .bs2_check_data(t1, t2)
  if (!.bs2_theta_valid(theta)) stop("Invalid theta.")
  f <- function(th) .bs2_loglik_theta(th, t1, t2)
  H <- .num_hessian(f, theta, h = .bs2_fd_steps(theta, rel_step = rel_step))
  J <- -H
  J <- (J + t(J)) / 2
  if (any(!is.finite(J))) stop("Observed information contains non-finite values.")
  diag(J) <- diag(J) + ridge
  J
}

# Gauss-Hermite quadrature for the MM covariance beta block.
.gh_cache <- new.env(parent = emptyenv())
.gh_stdnormal <- function(n = 15L) {
  key <- as.character(n)
  if (exists(key, envir = .gh_cache, inherits = FALSE)) return(get(key, envir = .gh_cache))
  J <- matrix(0, n, n)
  if (n > 1L) {
    off <- sqrt(seq_len(n - 1L) / 2)
    J[cbind(seq_len(n - 1L), 2:n)] <- off
    J[cbind(2:n, seq_len(n - 1L))] <- off
  }
  ev <- eigen(J, symmetric = TRUE)
  ord <- order(ev$values)
  nodes <- ev$values[ord]
  weights <- sqrt(pi) * ev$vectors[1, ord]^2
  out <- list(z = sqrt(2) * nodes, w = weights / sqrt(pi))
  assign(key, out, envir = .gh_cache)
  out
}

.bs_ratio <- function(z, alpha) {
  (alpha * z / 2 + sqrt((alpha * z / 2)^2 + 1))^2
}

psi_star_bs2 <- function(alpha, rho, gh_n = 15L) {
  rho <- .clip_rho(rho, eps = 1e-12)
  gh <- .gh_stdnormal(gh_n)
  u <- gh$z
  v <- gh$z
  w <- gh$w
  Z1 <- matrix(rep(u, times = length(v)), nrow = length(u), ncol = length(v))
  V <- matrix(rep(v, each = length(u)), nrow = length(u), ncol = length(v))
  Z2 <- rho * Z1 + sqrt(1 - rho^2) * V
  W <- outer(w, w)
  sum(W * .bs_ratio(Z1, alpha[1]) * .bs_ratio(Z2, alpha[2]))
}

mm_cov_bs2 <- function(alpha, beta, rho, gh_n = 15L) {
  alpha <- as.numeric(alpha); beta <- as.numeric(beta); rho <- .clip_rho(rho)
  a1 <- alpha[1]; a2 <- alpha[2]; b1 <- beta[1]; b2 <- beta[2]
  Xi <- matrix(0, 5, 5)

  Xi[1:2, 1:2] <- 0.5 * matrix(c(a1^2, a1 * a2 * rho^2,
                                   a1 * a2 * rho^2, a2^2), 2, 2)
  Xi[1:2, 5] <- rho * (1 - rho^2) / 2 * c(a1, a2)
  Xi[5, 1:2] <- Xi[1:2, 5]
  Xi[5, 5] <- (1 - rho^2)^2

  xi11 <- (a1 * b1)^2 * (3 * a1^2 + 4) / (a1^2 + 2)^2
  xi22 <- (a2 * b2)^2 * (3 * a2^2 + 4) / (a2^2 + 2)^2
  c_diff <- psi_star_bs2(alpha, rho, gh_n = gh_n) - psi_star_bs2(alpha, -rho, gh_n = gh_n)
  xi12 <- 2 * b1 * b2 * c_diff / ((a1^2 + 2) * (a2^2 + 2))
  Xi[3:4, 3:4] <- matrix(c(xi11, xi12, xi12, xi22), 2, 2)

  Xi
}
# Restricted profile maximum likelihood fits for H01-H04.
# Requires 00_helpers.R.

fit_bs2_restricted <- function(t1, t2, hypothesis, rho0 = 0, beta_start = NULL,
                               rho_eps = 1e-8,
                               optim_methods = c("BFGS", "Nelder-Mead"),
                               control = list(maxit = 1000)) {
  .bs2_check_data(t1, t2)
  r <- .restriction_matrix(hypothesis, rho0 = rho0)
  h <- r$hypothesis
  n <- length(t1)

  if (is.null(beta_start)) beta_start <- fit_bs2_mm(t1, t2)$beta
  beta_start <- as.numeric(beta_start)
  if (length(beta_start) == 1L) beta_start <- rep(beta_start, 2L)
  if (length(beta_start) != 2L || any(!is.finite(beta_start)) || any(beta_start <= 0)) {
    beta_start <- fit_bs2_mm(t1, t2)$beta
  }

  prof_from_eta <- function(eta) {
    if (h == "H01") {
      beta <- exp(eta)
      q1 <- a_bs_unit_alpha(t1, beta[1])
      q2 <- a_bs_unit_alpha(t2, beta[2])
      s11 <- sum(q1^2); s22 <- sum(q2^2); s12 <- sum(q1 * q2)
      alpha_scalar <- sqrt((s11 + s22) / (2 * n))
      rho <- if ((s11 + s22) > 0) 2 * s12 / (s11 + s22) else 0
      rho <- .clip_rho(rho, eps = rho_eps)
      alpha <- rep(alpha_scalar, 2L)
    } else if (h == "H02") {
      beta_scalar <- exp(eta[1])
      beta <- rep(beta_scalar, 2L)
      q1 <- a_bs_unit_alpha(t1, beta_scalar)
      q2 <- a_bs_unit_alpha(t2, beta_scalar)
      s11 <- sum(q1^2); s22 <- sum(q2^2); s12 <- sum(q1 * q2)
      if (s11 <= 0 || s22 <= 0) return(NULL)
      alpha <- c(sqrt(s11 / n), sqrt(s22 / n))
      rho <- s12 / sqrt(s11 * s22)
      rho <- .clip_rho(rho, eps = rho_eps)
    } else if (h == "H03") {
      beta_scalar <- exp(eta[1])
      beta <- rep(beta_scalar, 2L)
      q1 <- a_bs_unit_alpha(t1, beta_scalar)
      q2 <- a_bs_unit_alpha(t2, beta_scalar)
      s11 <- sum(q1^2); s22 <- sum(q2^2); s12 <- sum(q1 * q2)
      alpha_scalar <- sqrt((s11 + s22) / (2 * n))
      alpha <- rep(alpha_scalar, 2L)
      rho <- if ((s11 + s22) > 0) 2 * s12 / (s11 + s22) else 0
      rho <- .clip_rho(rho, eps = rho_eps)
    } else if (h == "H04") {
      beta <- exp(eta)
      rho <- .clip_rho(rho0, eps = rho_eps)
      q1 <- a_bs_unit_alpha(t1, beta[1])
      q2 <- a_bs_unit_alpha(t2, beta[2])
      s11 <- sum(q1^2); s22 <- sum(q2^2); s12 <- sum(q1 * q2)
      if (s11 <= 0 || s22 <= 0) return(NULL)
      # Exact profiled alpha for fixed rho and fixed beta.
      # For rho0 = 0 this reduces to alpha_j^2 = n^{-1} sum q_j^2.
      den <- n * (1 - rho^2)
      alpha1_sq <- (s11 - rho * s12 * sqrt(s11 / s22)) / den
      alpha2_sq <- (s22 - rho * s12 * sqrt(s22 / s11)) / den
      if (!is.finite(alpha1_sq) || !is.finite(alpha2_sq) || alpha1_sq <= 0 || alpha2_sq <= 0) return(NULL)
      alpha <- c(sqrt(alpha1_sq), sqrt(alpha2_sq))
    } else {
      stop("Unknown hypothesis.")
    }

    logLik <- ell_bs2(t1, t2, alpha = alpha, beta = beta, rho = rho)
    if (!is.finite(logLik)) return(NULL)
    list(alpha = alpha, beta = beta, rho = rho, logLik = logLik)
  }

  if (h %in% c("H01", "H04")) {
    start_eta <- log(beta_start)
  } else {
    start_eta <- log(exp(mean(log(beta_start))))
  }

  objective <- function(eta) {
    fit <- prof_from_eta(eta)
    if (is.null(fit) || !is.finite(fit$logLik)) return(.Machine$double.xmax / 100)
    -fit$logLik
  }

  opt <- .optim_best(start_eta, objective, methods = optim_methods, control = control)
  fit <- prof_from_eta(opt$par)
  if (is.null(fit)) stop("Restricted fit failed at the optimizer solution.")

  fit$hypothesis <- h
  fit$rho0 <- if (h == "H04") rho0 else NA_real_
  fit$convergence <- opt$convergence
  fit$message <- opt$message
  fit$objective <- opt$value
  fit$theta <- .fit_to_theta(fit)
  fit
}
# LR, score, Wald and MM-Wald tests for H01-H04.
# Requires 00_helpers.R and fit_bs2_restricted.R.

compute_bs2_tests <- function(t1, t2, fit_u = NULL, fit_r = NULL, fit_mm = NULL,
                              hypothesis, rho0 = 0, gh_n = 15L,
                              rel_step = 1e-5, return_matrices = FALSE) {
  .bs2_check_data(t1, t2)
  r <- .restriction_matrix(hypothesis, rho0 = rho0)
  h <- r$hypothesis
  A <- r$A
  q <- r$q
  df <- r$df
  n <- length(t1)

  if (is.null(fit_mm)) fit_mm <- fit_bs2_mm(t1, t2)
  if (is.null(fit_u)) fit_u <- fit_bs2_profile(t1, t2, beta_start = fit_mm$beta)
  if (is.null(fit_r)) fit_r <- fit_bs2_restricted(t1, t2, h, rho0 = rho0, beta_start = fit_u$beta)

  theta_u <- .fit_to_theta(fit_u)
  theta_r <- .fit_to_theta(fit_r)
  theta_m <- c(as.numeric(fit_mm$alpha), as.numeric(fit_mm$beta), as.numeric(fit_mm$rho))

  LR <- 2 * (fit_u$logLik - fit_r$logLik)
  if (is.finite(LR) && LR < 0 && abs(LR) < 1e-7) LR <- 0

  Wald <- NA_real_
  Score <- NA_real_
  MM_Wald <- NA_real_
  J_u <- J_r <- Cov_u <- Cov_r <- U_r <- NULL

  # Wald statistic based on observed information at the unrestricted MLE.
  wald_try <- try({
    J_u <- obs_info_bs2(t1, t2, theta_u, rel_step = rel_step)
    Cov_u <- .safe_solve(J_u)
    diff_u <- as.numeric(A %*% theta_u - q)
    V_u <- A %*% Cov_u %*% t(A)
    Wald <- as.numeric(t(diff_u) %*% .safe_solve(V_u) %*% diff_u)
    if (is.finite(Wald) && Wald < 0 && abs(Wald) < 1e-7) Wald <- 0
  }, silent = TRUE)

  # Score statistic based on observed information at the restricted MLE.
  score_try <- try({
    U_r <- score_bs2(t1, t2, theta_r, rel_step = rel_step)
    J_r <- obs_info_bs2(t1, t2, theta_r, rel_step = rel_step)
    Cov_r <- .safe_solve(J_r)
    Score <- as.numeric(t(U_r) %*% Cov_r %*% U_r)
    if (is.finite(Score) && Score < 0 && abs(Score) < 1e-7) Score <- 0
  }, silent = TRUE)

  # MM-Wald statistic using the asymptotic covariance matrix in Theorem 4.
  mm_try <- try({
    Xi_m <- mm_cov_bs2(alpha = fit_mm$alpha, beta = fit_mm$beta, rho = fit_mm$rho, gh_n = gh_n)
    diff_m <- as.numeric(A %*% theta_m - q)
    V_m <- A %*% Xi_m %*% t(A)
    MM_Wald <- as.numeric(n * t(diff_m) %*% .safe_solve(V_m) %*% diff_m)
    if (is.finite(MM_Wald) && MM_Wald < 0 && abs(MM_Wald) < 1e-7) MM_Wald <- 0
  }, silent = TRUE)

  out <- list(
    hypothesis = h,
    df = df,
    LR = LR,
    Score = Score,
    Wald = Wald,
    MM_Wald = MM_Wald,
    LR_p = pchisq(LR, df = df, lower.tail = FALSE),
    Score_p = pchisq(Score, df = df, lower.tail = FALSE),
    Wald_p = pchisq(Wald, df = df, lower.tail = FALSE),
    MM_Wald_p = pchisq(MM_Wald, df = df, lower.tail = FALSE),
    fit_u = fit_u,
    fit_r = fit_r,
    fit_mm = fit_mm,
    wald_ok = !inherits(wald_try, "try-error"),
    score_ok = !inherits(score_try, "try-error"),
    mm_wald_ok = !inherits(mm_try, "try-error")
  )

  if (return_matrices) {
    out$A <- A
    out$q <- q
    out$theta_u <- theta_u
    out$theta_r <- theta_r
    out$theta_m <- theta_m
    out$J_u <- J_u
    out$J_r <- J_r
    out$Cov_u <- Cov_u
    out$Cov_r <- Cov_r
    out$U_r <- U_r
  }

  out
}
# Parametric bootstrap calibration of the LR statistic under H01-H04.
# Requires 00_helpers.R, fit_bs2_restricted.R, and optionally compute_bs2_tests.R.

boot_lr_h0 <- function(t1, t2, hypothesis, B = 499L, rho0 = 0, seed = NULL,
                       beta_start = NULL, parallel = FALSE, cores = 1L,
                       verbose = FALSE) {
  .bs2_check_data(t1, t2)
  r <- .restriction_matrix(hypothesis, rho0 = rho0)
  h <- r$hypothesis
  n <- length(t1)
  if (!is.null(seed)) set.seed(seed)

  fit_mm <- fit_bs2_mm(t1, t2)
  if (is.null(beta_start)) beta_start <- fit_mm$beta
  fit_u <- fit_bs2_profile(t1, t2, beta_start = beta_start)
  fit_r <- fit_bs2_restricted(t1, t2, h, rho0 = rho0, beta_start = fit_u$beta)
  lr_obs <- 2 * (fit_u$logLik - fit_r$logLik)
  if (is.finite(lr_obs) && lr_obs < 0 && abs(lr_obs) < 1e-7) lr_obs <- 0

  one_boot <- function(b) {
    ans <- try({
      dat_b <- rbs2(n, alpha = fit_r$alpha, beta = fit_r$beta, rho = fit_r$rho)
      fit_mb <- fit_bs2_mm(dat_b[, 1], dat_b[, 2])
      fit_ub <- fit_bs2_profile(dat_b[, 1], dat_b[, 2], beta_start = fit_mb$beta)
      fit_rb <- fit_bs2_restricted(dat_b[, 1], dat_b[, 2], h, rho0 = rho0, beta_start = fit_ub$beta)
      lr <- 2 * (fit_ub$logLik - fit_rb$logLik)
      if (is.finite(lr) && lr < 0 && abs(lr) < 1e-7) lr <- 0
      lr
    }, silent = TRUE)
    if (inherits(ans, "try-error") || !is.finite(ans)) NA_real_ else as.numeric(ans)
  }

  if (isTRUE(parallel) && cores > 1L && .Platform$OS.type != "windows") {
    lr_star <- unlist(parallel::mclapply(seq_len(B), one_boot, mc.cores = cores, mc.set.seed = TRUE), use.names = FALSE)
  } else {
    lr_star <- vapply(seq_len(B), one_boot, numeric(1L))
  }

  lr_star_ok <- lr_star[is.finite(lr_star)]
  B_eff <- length(lr_star_ok)
  if (B_eff == 0L) stop("All bootstrap replications failed.")
  if (verbose && B_eff < B) message("Bootstrap failures: ", B - B_eff, " / ", B)

  list(
    hypothesis = h,
    df = r$df,
    lr_obs = lr_obs,
    lr_star = lr_star,
    B = B,
    B_eff = B_eff,
    failure_rate = 1 - B_eff / B,
    crit95 = as.numeric(quantile(lr_star_ok, 0.95, type = 1, names = FALSE)),
    p_boot = (1 + sum(lr_star_ok >= lr_obs)) / (B_eff + 1),
    p_asym = pchisq(lr_obs, df = r$df, lower.tail = FALSE),
    fit_u = fit_u,
    fit_r = fit_r
  )
}
# Cross-fitted conditional pseudo-likelihood for the BS2 dependence parameter.
# Requires 00_helpers.R.

fit_bs2_pl_crossfit <- function(t1, t2, K = 5L, folds = NULL, seed = NULL,
                                random_folds = TRUE, penalized = FALSE,
                                lambda_max = 50) {
  .bs2_check_data(t1, t2)
  n <- length(t1)
  K <- as.integer(K)
  if (K < 2L) stop("K must be at least 2.")
  K <- min(K, n)

  if (is.null(folds)) {
    idx <- seq_len(n)
    if (isTRUE(random_folds)) {
      if (!is.null(seed)) set.seed(seed)
      idx <- sample.int(n)
    }
    fold_id <- rep(seq_len(K), length.out = n)
    folds <- split(idx, fold_id)
  } else {
    if (length(folds) != K) K <- length(folds)
    folds <- lapply(folds, as.integer)
  }

  x_all <- numeric(0L)
  c_all <- numeric(0L)
  fold_fits <- vector("list", K)

  for (k in seq_len(K)) {
    test_idx <- folds[[k]]
    train_idx <- setdiff(seq_len(n), test_idx)
    if (length(train_idx) < 2L) stop("Each training fold must contain at least two observations.")

    fit_k <- fit_bs2_mm(t1[train_idx], t2[train_idx])
    fold_fits[[k]] <- fit_k

    y1 <- a_bs(t1[test_idx], fit_k$alpha[1], fit_k$beta[1])
    y2 <- a_bs(t2[test_idx], fit_k$alpha[2], fit_k$beta[2])
    x_all <- c(x_all, abs(y1 * y2))
    c_all <- c(c_all, as.numeric(y1 * y2 > 0))
  }

  fit <- .fit_conditional_logistic_slope(x_all, c_all, penalized = penalized, lambda_max = lambda_max)
  lambda_hat <- fit$lambda
  rho_hat <- fit$rho

  ell0 <- sum(c_all * 0 - log1pexp(0 * x_all))
  lr0 <- if (is.finite(lambda_hat)) 2 * (fit$logLik - ell0) else NA_real_

  pi_hat <- if (is.finite(lambda_hat)) plogis(lambda_hat * x_all) else rep(ifelse(lambda_hat > 0, 1, 0), length(x_all))
  A_hat <- mean(x_all^2 * pi_hat * (1 - pi_hat))
  B_hat <- mean(x_all^2 * (c_all - pi_hat)^2)
  se_lambda <- if (is.finite(A_hat) && A_hat > 0) sqrt((1 / length(x_all)) * B_hat / A_hat^2) else NA_real_
  gprime <- if (abs(rho_hat) < 1) (1 - rho_hat^2)^2 / (2 * (1 + rho_hat^2)) else NA_real_
  se_rho <- if (is.finite(gprime) && is.finite(se_lambda)) abs(gprime) * se_lambda else NA_real_

  score0 <- if (sum(x_all^2) > 0) (sum(x_all * (2 * c_all - 1))^2) / sum(x_all^2) else NA_real_

  list(lambda = lambda_hat, rho = rho_hat, se_lambda = se_lambda, se_rho = se_rho,
       logLik = fit$logLik, lr0 = lr0, score0 = score0, x = x_all, c = c_all,
       K = K, folds = folds, fold_fits = fold_fits,
       separation = fit$separation, boundary = fit$boundary, convergence = fit$convergence)
}
# Analytic score and expected Fisher information for LR/Score/Wald tests H01--H04.
# This block is integrated into the public core implementation; no external patch
# file is required.
# -----------------------------------------------------------------------------

J_bs <- function(alpha, gh_n = 25L) {
  alpha <- as.numeric(alpha)[1L]
  if (!is.finite(alpha) || alpha <= 0) stop("alpha must be positive.")
  gh <- .gh_stdnormal(gh_n)
  z <- gh$z
  w <- gh$w
  r <- .bs_ratio(z, alpha)
  sum(w * (r + 1)^(-2))
}

psi1_bs2 <- function(alpha, rho, gh_n = 25L) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 2L || any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop("alpha must be a positive vector of length 2.")
  }
  rho <- .clip_rho(rho, eps = 1e-12)
  gh <- .gh_stdnormal(gh_n)
  u <- gh$z
  v <- gh$z
  w <- gh$w

  Z1 <- matrix(rep(u, times = length(v)), nrow = length(u), ncol = length(v))
  V  <- matrix(rep(v, each  = length(u)), nrow = length(u), ncol = length(v))
  Z2 <- rho * Z1 + sqrt(1 - rho^2) * V
  W  <- outer(w, w)

  even_part <- sqrt(1 + (alpha[1]^2 * Z1^2) / 4) *
    sqrt(1 + (alpha[2]^2 * Z2^2) / 4)

  alpha[1] * alpha[2] * rho / 4 + sum(W * even_part)
}

fisher_info_bs2 <- function(alpha, beta, rho, gh_n = 25L, ridge = 0) {
  alpha <- as.numeric(alpha)
  beta  <- as.numeric(beta)
  if (length(alpha) != 2L || length(beta) != 2L) stop("alpha and beta must have length 2.")
  if (any(!is.finite(alpha)) || any(alpha <= 0)) stop("alpha must be positive.")
  if (any(!is.finite(beta)) || any(beta <= 0)) stop("beta must be positive.")
  rho <- .clip_rho(rho, eps = 1e-10)

  a1 <- alpha[1]; a2 <- alpha[2]
  b1 <- beta[1];  b2 <- beta[2]
  den <- 1 - rho^2

  I <- matrix(0, 5, 5)

  Dainv <- diag(1 / alpha, 2, 2)
  one <- matrix(1, 2, 1)
  Iaa <- (2 / den) * diag(1 / alpha^2, 2, 2) -
    (rho^2 / den) * (Dainv %*% one %*% t(one) %*% Dainv)

  Iar <- -(rho / den) * matrix(1 / alpha, ncol = 1)
  Irr <- (1 + rho^2) / den^2

  Jv <- c(J_bs(a1, gh_n = gh_n), J_bs(a2, gh_n = gh_n))
  Ibb <- matrix(0, 2, 2)
  Ibb[1, 1] <- ((rho^2 / 4 + 1 / a1^2) / den + Jv[1]) / b1^2
  Ibb[2, 2] <- ((rho^2 / 4 + 1 / a2^2) / den + Jv[2]) / b2^2

  psi_sum <- psi1_bs2(alpha, rho, gh_n = gh_n) + psi1_bs2(alpha, -rho, gh_n = gh_n)
  off <- -rho * psi_sum / (2 * den * a1 * a2 * b1 * b2)
  Ibb[1, 2] <- off
  Ibb[2, 1] <- off

  I[1:2, 1:2] <- Iaa
  I[1:2, 5]   <- Iar[, 1]
  I[5, 1:2]   <- Iar[, 1]
  I[5, 5]     <- Irr
  I[3:4, 3:4] <- Ibb

  I <- (I + t(I)) / 2
  if (ridge > 0) diag(I) <- diag(I) + ridge
  I
}

score_bs2_analytic <- function(t1, t2, theta) {
  .bs2_check_data(t1, t2)
  theta <- as.numeric(theta)
  if (!.bs2_theta_valid(theta)) stop("Invalid theta.")

  alpha <- theta[1:2]
  beta  <- theta[3:4]
  rho   <- .clip_rho(theta[5], eps = 1e-10)
  n <- length(t1)

  z1 <- a_bs(t1, alpha[1], beta[1])
  z2 <- a_bs(t2, alpha[2], beta[2])
  Z <- cbind(z1, z2)

  den <- 1 - rho^2
  Sinv <- matrix(c(1, -rho, -rho, 1), 2, 2) / den
  W <- Z %*% Sinv

  U_alpha <- c(
    (-n + sum(z1 * W[, 1])) / alpha[1],
    (-n + sum(z2 * W[, 2])) / alpha[2]
  )

  u1 <- sqrt(t1 / beta[1]); v1 <- sqrt(beta[1] / t1)
  u2 <- sqrt(t2 / beta[2]); v2 <- sqrt(beta[2] / t2)

  dlogA_db1 <- (beta[1] - t1) / (2 * beta[1] * (t1 + beta[1]))
  dlogA_db2 <- (beta[2] - t2) / (2 * beta[2] * (t2 + beta[2]))

  dz1_db1 <- -(u1 + v1) / (2 * alpha[1] * beta[1])
  dz2_db2 <- -(u2 + v2) / (2 * alpha[2] * beta[2])

  U_beta <- c(
    sum(dlogA_db1 - W[, 1] * dz1_db1),
    sum(dlogA_db2 - W[, 2] * dz2_db2)
  )

  Sigma_rho <- matrix(c(0, 1, 1, 0), 2, 2)
  Sinv_rho <- -Sinv %*% Sigma_rho %*% Sinv
  quad_rho <- rowSums((Z %*% Sinv_rho) * Z)
  U_rho <- n * rho / den - 0.5 * sum(quad_rho)

  c(U_alpha, U_beta, U_rho)
}

# Override the numerical finite-difference score with the analytic one.
score_bs2 <- function(t1, t2, theta, rel_step = NULL) {
  score_bs2_analytic(t1, t2, theta)
}

compute_bs2_tests <- function(t1, t2, fit_u = NULL, fit_r = NULL, fit_mm = NULL,
                              hypothesis, rho0 = 0, gh_n = 25L,
                              rel_step = 1e-5, return_matrices = FALSE,
                              score_type = c("efficient", "full")) {
  .bs2_check_data(t1, t2)
  score_type <- match.arg(score_type)
  r <- .restriction_matrix(hypothesis, rho0 = rho0)
  h <- r$hypothesis
  A <- r$A
  q <- r$q
  df <- r$df
  n <- length(t1)

  if (is.null(fit_mm)) fit_mm <- fit_bs2_mm(t1, t2)
  if (is.null(fit_u)) fit_u <- fit_bs2_profile(t1, t2, beta_start = fit_mm$beta)
  if (is.null(fit_r)) fit_r <- fit_bs2_restricted(t1, t2, h, rho0 = rho0, beta_start = fit_u$beta)

  theta_u <- .fit_to_theta(fit_u)
  theta_r <- .fit_to_theta(fit_r)
  theta_m <- c(as.numeric(fit_mm$alpha), as.numeric(fit_mm$beta), as.numeric(fit_mm$rho))

  LR <- 2 * (fit_u$logLik - fit_r$logLik)
  if (is.finite(LR) && LR < 0 && abs(LR) < 1e-7) LR <- 0

  Wald <- NA_real_
  Score <- NA_real_
  Score_full <- NA_real_
  MM_Wald <- NA_real_
  I_u <- I_r <- Omega_u <- Omega_r <- U_r <- NULL
  wald_ok <- score_ok <- mm_wald_ok <- FALSE

  wald_try <- try({
    I_u <- fisher_info_bs2(alpha = theta_u[1:2], beta = theta_u[3:4], rho = theta_u[5], gh_n = gh_n)
    Omega_u <- .safe_solve(I_u)
    diff_u <- as.numeric(A %*% theta_u - q)
    V_u <- A %*% Omega_u %*% t(A)
    Wald <- as.numeric(n * t(diff_u) %*% .safe_solve(V_u) %*% diff_u)
    if (is.finite(Wald) && Wald < 0 && abs(Wald) < 1e-7) Wald <- 0
    wald_ok <- TRUE
  }, silent = TRUE)

  score_try <- try({
    U_r <- score_bs2_analytic(t1, t2, theta_r)
    I_r <- fisher_info_bs2(alpha = theta_r[1:2], beta = theta_r[3:4], rho = theta_r[5], gh_n = gh_n)
    Omega_r <- .safe_solve(I_r)

    Score_full <- as.numeric(t(U_r) %*% Omega_r %*% U_r / n)

    AOAt <- A %*% Omega_r %*% t(A)
    P_eff <- Omega_r %*% t(A) %*% .safe_solve(AOAt) %*% A %*% Omega_r
    Score_eff <- as.numeric(t(U_r) %*% P_eff %*% U_r / n)

    Score <- if (score_type == "efficient") Score_eff else Score_full
    if (is.finite(Score) && Score < 0 && abs(Score) < 1e-7) Score <- 0
    if (is.finite(Score_full) && Score_full < 0 && abs(Score_full) < 1e-7) Score_full <- 0
    score_ok <- TRUE
  }, silent = TRUE)

  mm_try <- try({
    Xi_m <- mm_cov_bs2(alpha = fit_mm$alpha, beta = fit_mm$beta, rho = fit_mm$rho, gh_n = gh_n)
    diff_m <- as.numeric(A %*% theta_m - q)
    V_m <- A %*% Xi_m %*% t(A)
    MM_Wald <- as.numeric(n * t(diff_m) %*% .safe_solve(V_m) %*% diff_m)
    if (is.finite(MM_Wald) && MM_Wald < 0 && abs(MM_Wald) < 1e-7) MM_Wald <- 0
    mm_wald_ok <- TRUE
  }, silent = TRUE)

  out <- list(
    hypothesis = h,
    df = df,
    LR = LR,
    Score = Score,
    Score_full = Score_full,
    Wald = Wald,
    MM_Wald = MM_Wald,
    LR_p = pchisq(LR, df = df, lower.tail = FALSE),
    Score_p = pchisq(Score, df = df, lower.tail = FALSE),
    Score_full_p = pchisq(Score_full, df = df, lower.tail = FALSE),
    Wald_p = pchisq(Wald, df = df, lower.tail = FALSE),
    MM_Wald_p = pchisq(MM_Wald, df = df, lower.tail = FALSE),
    fit_u = fit_u,
    fit_r = fit_r,
    fit_mm = fit_mm,
    wald_ok = wald_ok && !inherits(wald_try, "try-error"),
    score_ok = score_ok && !inherits(score_try, "try-error"),
    mm_wald_ok = mm_wald_ok && !inherits(mm_try, "try-error"),
    score_type = score_type
  )

  if (return_matrices) {
    out$A <- A
    out$q <- q
    out$theta_u <- theta_u
    out$theta_r <- theta_r
    out$theta_m <- theta_m
    out$I_u <- I_u
    out$I_r <- I_r
    out$Omega_u <- Omega_u
    out$Omega_r <- Omega_r
    out$U_r <- U_r
  }

  out
}
