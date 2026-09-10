#!/usr/bin/env Rscript
#
# Run the toolkit test suite:  Rscript tests/run_tests.R
# Loads the toolkit from the directory above this file, then runs
# tests/testthat with testthat and stops with a non-zero status on failure.

script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
home <- if (length(script) > 0) normalizePath(file.path(dirname(script), "..")) else normalizePath(".")
options(wca.quiet = TRUE)
source(file.path(home, "load.R"))
suppressPackageStartupMessages(library(testthat))
res <- testthat::test_dir(file.path(home, "tests", "testthat"), reporter = "summary",
  stop_on_failure = FALSE, stop_on_warning = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n%d tests: %d passed, %d failed, %d skipped\n",
  nrow(df), sum(df$passed), sum(df$failed), sum(df$skipped)))
if (sum(df$failed) > 0 || any(df$error)) quit(status = 1)
