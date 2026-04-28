#!/usr/bin/env Rscript
# Install the R packages required by the BS2 reproducibility bundle.

repos <- getOption("repos")
if (is.null(repos) || identical(repos[["CRAN"]], "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

required <- c(
  "data.table",
  "future",
  "future.apply",
  "parallelly",
  "RhpcBLASctl",
  "ggplot2"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing)
} else {
  message("All required packages are already installed.")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Optional package 'renv' is not installed. Install it if you want a lockfile:")
  message("  install.packages('renv'); renv::init(bare = TRUE); renv::snapshot()")
}

message("\nInstalled package versions:")
for (pkg in required) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  %-15s %s", pkg, as.character(utils::packageVersion(pkg))))
  }
}
