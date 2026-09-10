#!/usr/bin/env Rscript
#
# wcaUpdateManifests.R: add new pairs to, or update reruns in, the SNV and
# SV Tempo manifests from a list of Tempo output roots.
#
#   Rscript $WCA_HOME/bin/wcaUpdateManifests.R \
#     <snv_manifest.csv> <sv_manifest.csv> <roots_file> [out_dir] \
#     [--stamp=yymmdd] [--force] [--dedup[=first|newest|oldest]] \
#     [--prefer=<regex>]
#
# roots_file: one Tempo output root per line (a directory containing
# out/<run>/somatic/<PAIR>/...); blank lines and # comments are ignored;
# an optional ",ProjNo" after the path overrides the project number that
# is otherwise taken from the run directory (Proj_ prefix removed).
#
# Writes tempoSNVManifest_<stamp>.csv, tempoSVManifest_<stamp>.csv,
# manifestUpdate_<stamp>.md and manifestUpdate_<stamp>.csv (one row per
# added or updated tumor, old and new values side by side) to out_dir
# (default: the directory of the SNV manifest). Refuses to overwrite
# existing outputs unless --force.
#
# Tumors listed more than once are always reported, in the run log and in
# the report. --dedup drops the rows that name a file already listed for
# the same tumor; every dropped row is listed in the report and the csv. A
# tumor listed with two different files is never dropped, whatever --dedup
# says, because only a person can decide which run is current.
#
# Of a redundant group the copy under out/ is kept, since that is where
# Tempo writes and a second path for the same md5 is a copy of it. Change
# that with --prefer=<regex on PATH>, or disable it with --prefer= and
# choose by position or age instead: --dedup keeps the first row of the
# group, --dedup=newest / --dedup=oldest the newest or oldest file.

argv <- commandArgs(trailingOnly = TRUE)
flags <- argv[startsWith(argv, "--")]
args <- argv[!startsWith(argv, "--")]
if (length(args) < 3) {
  stop("usage: wcaUpdateManifests.R <snv_manifest> <sv_manifest> <roots_file> [out_dir] [--stamp=yymmdd] [--force] [--dedup[=first|newest|oldest]]")
}
flag_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), flags, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[[1]])
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
options(wca.quiet = TRUE)
source(file.path(dirname(script), "..", "load.R"))
options(wca.quiet = FALSE)

inputs <- c(snv = args[[1]], sv = args[[2]])
roots_file <- args[[3]]
out_dir <- if (length(args) >= 4) args[[4]] else as.character(fs::path_dir(inputs[["snv"]]))
inputs <- set_names(as.character(fs::path_abs(inputs)), names(inputs))
out_dir <- as.character(fs::path_abs(out_dir))
stamp <- flag_value("stamp", date_stamp())
force <- "--force" %in% flags
dedup <- any(flags == "--dedup" | startsWith(flags, "--dedup="))
keep <- flag_value("dedup", "first")
if (!keep %in% c("first", "newest", "oldest")) {
  wca_abort("--dedup must be one of first, newest, oldest (got {keep})")
}
prefer <- flag_value("prefer", "/out/")

for (f in c(inputs, roots_file)) if (!fs::file_exists(f)) wca_abort("Input not found: {f}")

old <- list(snv = read_manifest(inputs[["snv"]]), sv = read_manifest(inputs[["sv"]]))
wca_msg("Read {nrow(old$snv)} SNV and {nrow(old$sv)} SV manifest rows")

root_spec <- read_list_file(roots_file) |>
  tibble(line = _) |>
  tidyr::separate_wider_delim(line, ",", names = c("root", "proj_no"), too_few = "align_start") |>
  mutate(across(everything(), str_trim), proj_no = na_if(proj_no, ""))
missing_roots <- root_spec$root[!fs::dir_exists(root_spec$root)]
if (length(missing_roots) > 0) wca_abort("Roots not found: {str_c(missing_roots, collapse = ', ')}")

scans <- purrr::pmap(root_spec, \(root, proj_no) {
  sc <- scan_tempo_outputs(root, proj_no = if (is.na(proj_no)) NULL else proj_no)
  wca_msg("  {root}: {nrow(sc$snv)} MAF, {nrow(sc$sv)} bedpe")
  sc
})
new <- list(
  snv = purrr::map(scans, "snv") |> bind_rows(),
  sv = purrr::map(scans, "sv") |> bind_rows()
)

updates <- list(
  snv = update_manifest(old$snv, new$snv, dedup = dedup, keep = keep, prefer = prefer),
  sv = update_manifest(old$sv, new$sv, dedup = dedup, keep = keep, prefer = prefer)
)
keep_rule <- if (nzchar(prefer)) {
  str_c("the path matching ", prefer, ", else the ", keep)
} else {
  str_c("the ", keep)
}
brief <- function(x, n = 8) {
  if (length(x) <= n) {
    str_c(x, collapse = ", ")
  } else {
    str_c(str_c(head(x, n), collapse = ", "), ", ... (", length(x), " total)")
  }
}
for (nm in names(updates)) {
  ch <- updates[[nm]]$changes
  wca_msg("  {nm}: {sum(ch$status == 'added')} added, {sum(ch$status == 'updated')} updated, ",
    "{sum(ch$status == 'unchanged')} unchanged, {sum(ch$status == 'kept')} kept")
  removed <- sum(ch$status == "removed")
  if (removed > 0) wca_msg("  {nm}: dropped {removed} duplicate row(s), keeping {keep_rule}")
  m <- updates[[nm]]$manifest
  left <- updates[[nm]]$duplicates
  exact <- left |> filter(dup == "exact") |> distinct(TID) |> pull(TID)
  conflict <- left |> filter(dup == "conflict") |> distinct(TID) |> pull(TID)
  if (length(exact) > 0) {
    wca_warn("  {nm}: {length(exact)} TID listed more than once with the same file ",
      "(rerun with --dedup to drop the extra rows): {brief(exact)}")
  }
  if (length(conflict) > 0) {
    wca_warn("  {nm}: {length(conflict)} TID listed with different files (resolve by hand): {brief(conflict)}")
  }
  gone <- m$PATH[!fs::file_exists(m$PATH)]
  if (length(gone) > 0) {
    wca_warn("  {nm}: {length(gone)} PATH entries do not exist (first: {gone[[1]]})")
  }
}

outputs <- write_manifests(list(snv = updates$snv$manifest, sv = updates$sv$manifest),
  out_dir, stamp = stamp, overwrite = force)
report <- fs::path(out_dir, str_glue("manifestUpdate_{stamp}.md"))
if (fs::file_exists(report) && !force) wca_abort("Report exists (use --force): {report}")
readr::write_lines(
  manifest_update_report(updates, inputs = inputs, roots = root_spec$root, outputs = outputs),
  report)
outputs[["report"]] <- as.character(report)

changes <- fs::path(out_dir, str_glue("manifestUpdate_{stamp}.csv"))
if (fs::file_exists(changes) && !force) wca_abort("Change table exists (use --force): {changes}")
readr::write_csv(manifest_update_changes(updates), changes)
outputs[["changes"]] <- as.character(changes)

for (nm in names(outputs)) wca_msg("wrote {outputs[[nm]]}")
