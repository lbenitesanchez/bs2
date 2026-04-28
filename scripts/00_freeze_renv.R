#!/usr/bin/env Rscript
# Create or update renv.lock after package installation.
# Run this on the Ubuntu machine used to generate final manuscript results.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
}
renv::snapshot(prompt = FALSE)
cat("Updated renv.lock\n")
