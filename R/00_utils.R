# 00_utils.R
#
# Small helpers shared by every other file. Nothing here touches the file
# system on load, and nothing depends on ~/.Rprofile.

#' Toolkit home directory
#'
#' @return The directory that holds `load.R`, taken from option `wca.home`
#'   (set by `load.R`) or the `WCA_HOME` environment variable.
wca_home <- function() {
  home <- getOption("wca.home", Sys.getenv("WCA_HOME", unset = ""))
  if (!nzchar(home)) {
    wca_abort("Toolkit home is not set. Source load.R or set WCA_HOME.")
  }
  home
}

#' Path to a file inside the toolkit
#'
#' @param ... Path components below the toolkit home.
#' @return A single path.
wca_file <- function(...) {
  fs::path(wca_home(), ...)
}

#' Toolkit version string
#'
#' @param with_sha Append the short git SHA when the toolkit is a git checkout.
#' @return `"0.1.0"` or `"0.1.0 (abc1234)"`.
toolkit_version <- function(with_sha = TRUE) {
  version <- readr::read_lines(wca_file("VERSION"), n_max = 1)
  sha <- if (with_sha) toolkit_git_sha() else NA_character_
  if (is.na(sha)) version else str_glue("{version} ({sha})")
}

#' Short git SHA of the toolkit checkout
#'
#' @return The SHA, or `NA` when git or the repository is unavailable.
toolkit_git_sha <- function() {
  home <- wca_home()
  if (!fs::dir_exists(fs::path(home, ".git"))) return(NA_character_)
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(home), "rev-parse", "--short", "HEAD"),
      stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  ))
  if (length(out) == 0 || !is.null(attr(out, "status"))) NA_character_ else out[[1]]
}

#' Date stamp for output file names
#'
#' @param date A `Date`; defaults to today.
#' @return A `yymmdd` string such as `"260909"`.
date_stamp <- function(date = Sys.Date()) {
  format(date, "%y%m%d")
}

#' Convert character columns to their natural types without messages
#'
#' @param x A tibble read with `col_types = cols(.default = "c")`.
#' @return The same tibble with columns type-converted.
type_convert_silent <- function(x) {
  suppressMessages(readr::type_convert(x, guess_integer = TRUE))
}

#' Stop with a toolkit-prefixed message
#'
#' @param ... Passed to `stringr::str_glue()`, evaluated in the caller.
#' @param .envir Environment for glue interpolation.
wca_abort <- function(..., .envir = parent.frame()) {
  rlang::abort(as.character(str_glue(..., .envir = .envir)), class = "wca_error")
}

#' Warn with a toolkit-prefixed message
#'
#' @inheritParams wca_abort
wca_warn <- function(..., .envir = parent.frame()) {
  rlang::warn(as.character(str_glue(..., .envir = .envir)), class = "wca_warning")
}

#' Print a progress line
#'
#' Output goes to stdout so it lands in run logs; suppress with
#' `options(wca.quiet = TRUE)`.
#'
#' @inheritParams wca_abort
wca_msg <- function(..., .envir = parent.frame()) {
  if (isTRUE(getOption("wca.quiet", FALSE))) return(invisible())
  cat(as.character(str_glue(..., .envir = .envir)), "\n", sep = "")
  invisible()
}

#' Assert that a data frame has the given columns
#'
#' @param x A data frame.
#' @param cols Required column names.
#' @param what Label used in the error message.
require_cols <- function(x, cols, what = "table") {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    wca_abort("{what} is missing columns: {str_c(missing, collapse = ', ')}")
  }
  invisible(x)
}

#' Read a whitespace-delimited list file into a character vector
#'
#' @param file Path to a file with one token per line; blank lines and lines
#'   starting with `#` are ignored.
#' @return Character vector.
read_list_file <- function(file) {
  readr::read_lines(file) |>
    str_trim() |>
    purrr::discard(\(x) x == "" | str_starts(x, "#"))
}

#' Null-coalescing helper
#'
#' @param a,b Values; returns `b` when `a` is `NULL`.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Split a Tempo pair name into tumor and normal ids
#'
#' @param pair Character vector of `TUMOR__NORMAL` names.
#' @return Tibble with columns `Sample` and `NID`.
split_pair_name <- function(pair) {
  tibble(pair = pair) |>
    tidyr::separate_wider_delim(
      pair, delim = "__", names = c("Sample", "NID"),
      too_many = "merge", too_few = "align_start"
    )
}
