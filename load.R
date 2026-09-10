# load.R
#
# Source this file to load the WGSCohortAnalysis toolkit into the current
# session:
#
#   source(file.path(Sys.getenv("WCA_HOME"), "load.R"))
#
# It attaches the tidyverse packages the toolkit uses, sources every file in
# R/ in name order, records the toolkit home in option `wca.home`, and prints
# the version and git SHA. Set `options(wca.quiet = TRUE)` before sourcing to
# suppress the banner.

wca_load <- function(home = NULL, quiet = getOption("wca.quiet", FALSE)) {
  if (is.null(home)) {
    home <- Sys.getenv("WCA_HOME", unset = "")
    if (!nzchar(home)) {
      stop("wca_load(): cannot locate the toolkit; pass home = or set WCA_HOME")
    }
  }
  home <- normalizePath(home, mustWork = TRUE)

  deps <- c(
    "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble", "forcats",
    "ggplot2", "fs", "glue", "vroom", "rlang"
  )
  missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("wca_load(): missing packages: ", paste(missing, collapse = ", "))
  }
  suppressPackageStartupMessages(
    invisible(lapply(deps, library, character.only = TRUE, warn.conflicts = FALSE))
  )

  options(wca.home = home)

  r_files <- sort(list.files(file.path(home, "R"), pattern = "\\.R$", full.names = TRUE))
  for (f in r_files) sys.source(f, envir = globalenv())

  if (!quiet) {
    cat(sprintf("WGSCohortAnalysis %s loaded from %s\n", toolkit_version(), home))
  }
  invisible(home)
}

local({
  this_file <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) {
      this_file <- f
      break
    }
  }
  home <- if (!is.null(this_file)) dirname(normalizePath(this_file)) else NULL
  wca_load(home = home)
})
