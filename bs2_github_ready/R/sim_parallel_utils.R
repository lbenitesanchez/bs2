# R/sim_parallel_utils.R
# -----------------------------------------------------------------------------
# Parallel execution utilities for the BS2 Monte Carlo simulations.
# Source this file after installing: data.table, future, future.apply, parallelly.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it with install.packages('data.table').")
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required. Install it with install.packages('future').")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required. Install it with install.packages('future.apply').")
  }
  if (!requireNamespace("parallelly", quietly = TRUE)) {
    stop("Package 'parallelly' is required. Install it with install.packages('parallelly').")
  }
  library(data.table)
  library(future)
  library(future.apply)
  library(parallelly)
})

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(pat, "", hit[[1L]])
}

get_int_arg <- function(name, default) {
  as.integer(get_arg(name, as.character(default)))
}

get_num_arg <- function(name, default) {
  as.numeric(get_arg(name, as.character(default)))
}

get_bool_arg <- function(name, default = FALSE) {
  val <- tolower(as.character(get_arg(name, as.character(default))))
  val %in% c("true", "t", "yes", "y", "1")
}

pin_numeric_threads <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }

  invisible(TRUE)
}

setup_parallel <- function(cores, backend = "multicore") {
  pin_numeric_threads()

  cores <- max(1L, as.integer(cores))

  if (backend == "multicore" && parallelly::supportsMulticore()) {
    future::plan(future::multicore, workers = cores)
  } else if (backend == "multisession") {
    future::plan(future::multisession, workers = cores)
  } else if (backend == "sequential" || cores == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = cores)
  }

  RNGkind("L'Ecuyer-CMRG")
  invisible(TRUE)
}

reset_parallel <- function() {
  future::plan(future::sequential)
  invisible(TRUE)
}

task_seed <- function(base_seed, cell_id, chunk_id) {
  s <- as.numeric(base_seed) +
    1000003 * as.numeric(cell_id) +
    9176 * as.numeric(chunk_id)
  s <- s %% .Machine$integer.max
  as.integer(s + 1L)
}

make_rep_tasks <- function(design, R, chunk_size, out_dir, prefix) {
  design <- as.data.table(design)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ans <- lapply(seq_len(nrow(design)), function(i) {
    starts <- seq.int(1L, R, by = chunk_size)
    ends <- pmin(starts + chunk_size - 1L, R)
    design_cell_id <- if ("cell_id" %in% names(design)) design$cell_id[i] else i

    data.table(
      design_row = i,
      cell_id = as.integer(design_cell_id),
      chunk_id = seq_along(starts),
      rep_start = starts,
      rep_end = ends,
      out_file = file.path(
        out_dir,
        sprintf("%s_cell%04d_chunk%04d.rds",
                prefix, as.integer(design_cell_id), seq_along(starts))
      )
    )
  })

  rbindlist(ans)
}

safe_eval <- function(expr) {
  tryCatch(expr, error = function(e) {
    structure(list(error = conditionMessage(e)), class = "try-error")
  })
}

fit_ok <- function(x) {
  !is.null(x) && !inherits(x, "try-error")
}

num_or_na <- function(x, name, idx = NULL) {
  if (is.null(x) || inherits(x, "try-error")) return(NA_real_)
  if (is.null(x[[name]])) return(NA_real_)
  val <- x[[name]]
  if (!is.null(idx)) {
    if (length(val) < idx) return(NA_real_)
    val <- val[idx]
  }
  suppressWarnings(as.numeric(val)[1L])
}

bool_or_na <- function(x, name) {
  if (is.null(x) || inherits(x, "try-error")) return(NA)
  if (is.null(x[[name]])) return(NA)
  as.logical(x[[name]][1L])
}

message_header <- function(...) {
  msg <- paste0(...)
  message(strrep("=", 72))
  message(msg)
  message(strrep("=", 72))
}

ensure_dirs <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

read_rds_dir <- function(path, pattern = "\\.rds$") {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) return(data.table())
  rbindlist(lapply(files, readRDS), fill = TRUE, use.names = TRUE)
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
