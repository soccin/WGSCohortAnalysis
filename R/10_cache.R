# 10_cache.R
#
# Content-addressed read cache. A parsed file is stored as an rds keyed by
# the md5 of the source file, the reader name and version, and a digest of
# the reader arguments, so editing the source, the reader, or the filter
# invalidates the entry. The cache lives at the project root
# (`cache/reads/`) outside the run axis because parsed reads are shared by
# every run of a project.

#' Bump this when a reader's output changes shape
#'
#' 2: readers no longer treat "," as a number grouping mark, so columns of
#' comma-separated pairs such as manta PR/SR and delly CIPOS/CIEND come back
#' as text instead of a run-together number. Entries written by version 1 are
#' unreachable after this bump, which is the point.
WCA_READER_VERSION <- "2"

#' Read a file through the cache
#'
#' @param file Source path.
#' @param reader Function `reader(file, ...)`.
#' @param ... Reader arguments; functions among them are digested by their
#'   `wca_filter` attribute when present, else by their deparsed body.
#' @param cache_dir Cache directory; created if needed.
#' @param reader_name Name used in the key; defaults to the deparsed
#'   `reader` argument.
#' @param refresh Ignore an existing entry.
#' @return The reader's result.
cached_read <- function(file, reader, ..., cache_dir, reader_name = NULL, refresh = FALSE) {
  reader_name <- reader_name %||% deparse1(substitute(reader))
  fs::dir_create(cache_dir)
  args <- list(...)
  key <- cache_key(file, reader_name, args)
  path <- fs::path(cache_dir, str_c(key, ".rds"))
  if (!refresh && fs::file_exists(path)) {
    return(readRDS(path))
  }
  out <- reader(file, ...)
  saveRDS(out, path, compress = "xz")
  out
}

#' Cache key for a file, reader and argument list
#'
#' @inheritParams cached_read
#' @param args List of reader arguments.
#' @return A 32-character hex string.
cache_key <- function(file, reader_name, args = list()) {
  file_md5 <- unname(tools::md5sum(file))
  norm_args <- purrr::map(args, \(a) {
    if (is.function(a)) attr(a, "wca_filter") %||% deparse1(body(a)) else a
  })
  digest::digest(list(file_md5, reader_name, WCA_READER_VERSION, norm_args), algo = "md5")
}

#' Remove every entry from a cache directory
#'
#' @param cache_dir Cache directory.
#' @return Number of files removed, invisibly.
cache_clear <- function(cache_dir) {
  if (!fs::dir_exists(cache_dir)) return(invisible(0L))
  files <- fs::dir_ls(cache_dir, glob = "*.rds")
  fs::file_delete(files)
  invisible(length(files))
}
