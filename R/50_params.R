# 50_params.R
#
# Project parameters and the run axis. An analysis project follows
# https://github.com/soccin/r-proj-org: `00.PARAMS.yml` at the root,
# stage outputs under `cache/<run>/NN_stage/`, deliverables under
# `results/<run>/`. The read cache (`cache/reads/`) sits outside the run
# axis on purpose: parsed Tempo files are shared by every run.

#' Default parameters
#'
#' @return Nested list mirroring `templates/project/00.PARAMS.yml`.
wca_default_params <- function() {
  list(
    run = "run01",
    toolkit = NULL,
    prefix = "wca",
    genome = "hg19",
    manifests = list(snv = NULL, sv = NULL),
    cohort = list(projects_file = NULL, include_ucode = NULL, include_projects = NULL,
      exclude_tid_regex = NULL, require_both = TRUE),
    snv = list(keep_5p_flank = "TERT", min_t_alt_count = 0, min_t_depth = 0, min_vaf = 0),
    cnv = list(max_abs_diplogr = 1.5, require_facets_qc = FALSE, autosomes_only = TRUE,
      drop_diploid = TRUE, diploid_rule = "state", call_preset = "facets", arm_min_frac = 0.5),
    sv = list(min_callers_pass = 0),
    report = list(pct_cutoff = list(snv = 0.10, sv = 0.10, cnv = 0.40),
      genes_of_interest = character(), top_n_oncoprint = 40)
  )
}

#' Read a yaml file where only true and false are booleans
#'
#' The yaml package follows YAML 1.1, where `y`, `n`, `yes`, `no`, `on` and
#' `off`, in any case, are booleans. In this toolkit's parameters those are
#' data: `Y` is tyrosine and a chromosome, `N` is asparagine and an unknown
#' base, so `ref_aa: Y` would arrive as `TRUE`. Here, as in YAML 1.2, only
#' `true` and `false` (any case) are logical; every other such word comes
#' back as the text written.
#'
#' @param file Path to the yaml.
#' @return The parsed yaml.
wca_read_yaml <- function(file) {
  text_unless_true_false <- function(x) switch(tolower(x), true = TRUE, false = FALSE, x)
  yaml::read_yaml(file, handlers = list(
    "bool#yes" = text_unless_true_false,
    "bool#no" = text_unless_true_false
  ))
}

#' Stop when a logical parameter is not true or false
#'
#' Under `wca_read_yaml()` a parameter written `yes` or `on` is text, not
#' `TRUE`. Code testing it with `isTRUE()` would read that as false without a
#' word, so every key whose default is logical must come back logical.
#'
#' @param p Parameters after merging with the defaults.
#' @param defaults The matching level of `wca_default_params()`.
#' @param path Keys above this level, for the message.
#' @return `p`, invisibly; stops with a `wca_user_error` otherwise.
check_logical_params <- function(p, defaults = wca_default_params(), path = character()) {
  purrr::iwalk(defaults, \(d, key) {
    v <- p[[key]]
    if (is.list(d) && is.list(v)) {
      check_logical_params(v, d, c(path, key))
    } else if (is.logical(d) && !is.null(v) && !(is.logical(v) && length(v) == 1 && !is.na(v))) {
      name <- str_c(c(path, key), collapse = ".")
      wca_stop_user(
        "{name} must be true or false",
        detail = c(
          str_glue("00.PARAMS.yml has {name}: {str_c(format(v), collapse = ', ')}"),
          "Words such as yes, no, on, off, y and n are read as text here, because Y and N are amino acids and bases."
        ),
        fix = str_glue("write {name}: true or {name}: false")
      )
    }
  })
  invisible(p)
}

#' Read a project's 00.PARAMS.yml
#'
#' Missing keys take their defaults; `toolkit` falls back to `WCA_HOME` and
#' then to the loaded toolkit. Relative paths are kept relative to the
#' project root (`root`). The yaml is read by `wca_read_yaml()`, so only
#' `true` and `false` are booleans, and a logical parameter written any
#' other way stops the read (`check_logical_params()`).
#'
#' @param file Path to the yaml.
#' @param root Project root; defaults to the yaml's directory.
#' @return Nested list with `root` added.
wca_read_params <- function(file = "00.PARAMS.yml", root = NULL) {
  if (!fs::file_exists(file)) wca_abort("Parameters file not found: {file}")
  root <- root %||% as.character(fs::path_dir(fs::path_abs(file)))
  raw <- wca_read_yaml(file)
  p <- modifyList(wca_default_params(), raw, keep.null = TRUE)
  check_logical_params(p)
  p$root <- root
  if (is.null(p$toolkit) || !nzchar(p$toolkit)) {
    p$toolkit <- Sys.getenv("WCA_HOME", unset = getOption("wca.home", ""))
  }
  if (!is.null(p$cohort$include_ucode)) p$cohort$include_ucode <- as.character(p$cohort$include_ucode)
  if (!is.null(p$cohort$include_projects)) p$cohort$include_projects <- as.character(p$cohort$include_projects)
  p$report$genes_of_interest <- as.character(p$report$genes_of_interest %||% character())
  p$snv$keep_5p_flank <- as.character(p$snv$keep_5p_flank %||% character())
  p
}

#' Resolve a project-relative path
#'
#' @param params From `wca_read_params()`.
#' @param ... Path components.
#' @return Absolute path.
project_path <- function(params, ...) {
  as.character(fs::path_abs(fs::path(...), start = params$root))
}

#' Directory for a stage's cache outputs (`cache/<run>/<stage>/`)
#'
#' @param params From `wca_read_params()`.
#' @param stage Stage name such as `"02_read"`.
#' @param create Create the directory.
#' @return Path.
cache_dir <- function(params, stage, create = TRUE) {
  d <- project_path(params, "cache", params$run, stage)
  if (create) fs::dir_create(d)
  d
}

#' Directory for the shared read cache (`cache/reads/`)
#'
#' @inheritParams cache_dir
#' @return Path.
reads_cache_dir <- function(params, create = TRUE) {
  d <- project_path(params, "cache", "reads")
  if (create) fs::dir_create(d)
  d
}

#' Directory under `results/<run>/`
#'
#' @inheritParams cache_dir
#' @param ... Sub-directories.
#' @return Path.
results_dir <- function(params, ..., create = TRUE) {
  d <- project_path(params, "results", params$run, ...)
  if (create) fs::dir_create(d)
  d
}

#' Save a stage object as rds and log it
#'
#' @param x Object.
#' @param params From `wca_read_params()`.
#' @param stage Stage name.
#' @param name File stem.
#' @return The path, invisibly.
stage_save <- function(x, params, stage, name) {
  path <- fs::path(cache_dir(params, stage), str_c(name, ".rds"))
  saveRDS(x, path)
  n <- if (is.data.frame(x)) str_glue("{nrow(x)} rows") else class(x)[[1]]
  wca_msg("  saved {stage}/{name}.rds ({n})")
  invisible(path)
}

#' Load a stage object saved by `stage_save()`
#'
#' @inheritParams stage_save
#' @return The object.
stage_load <- function(params, stage, name) {
  path <- fs::path(cache_dir(params, stage, create = FALSE), str_c(name, ".rds"))
  if (!fs::file_exists(path)) wca_abort("Missing stage output {path}; run the earlier stage first")
  readRDS(path)
}

#' Print a stage banner
#'
#' @param stage Stage name.
#' @param params From `wca_read_params()`.
stage_banner <- function(stage, params) {
  wca_msg("\n== {stage} [{params$run}] {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} ==")
}

#' Write TOOLKIT_VERSION and copy 00.PARAMS.yml into results/<run>/
#'
#' @param params From `wca_read_params()`.
#' @return Paths written, invisibly.
record_run_provenance <- function(params) {
  out <- results_dir(params)
  version_file <- fs::path(out, "TOOLKIT_VERSION")
  readr::write_lines(c(
    str_glue("WGSCohortAnalysis {toolkit_version()}"),
    str_glue("home: {wca_home()}"),
    str_glue("run: {params$run}"),
    str_glue("date: {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"),
    str_glue("R: {R.version.string}")
  ), version_file)
  params_src <- project_path(params, "00.PARAMS.yml")
  params_dst <- fs::path(out, "00.PARAMS.yml")
  if (fs::file_exists(params_src)) fs::file_copy(params_src, params_dst, overwrite = TRUE)
  invisible(c(version_file, params_dst))
}
