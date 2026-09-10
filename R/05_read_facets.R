# 05_read_facets.R
#
# FACETS readers. FACETS codes chromosome X as 23 and never emits Y, so every
# reader passes `chrom` through `normalize_chrom()`; an autosome filter must
# therefore be written against "1".."22" (see `autosomes()`).

facets_tsv <- function(file) {
  read_tsv(file, col_types = cols(.default = "c"), na = c("", "NA"), progress = FALSE) |>
    type_convert_silent()
}

#' Read a FACETS gene-level table
#'
#' @param file `*.gene_level.txt`.
#' @param drop_diploid Drop rows whose `cn_state` is `DIPLOID` (kept rows
#'   still include `tcn == 2` states such as CNLOH).
#' @return Tibble with `chrom` normalized and `sample` split into
#'   `Sample`/`NID`.
read_facets_gene <- function(file, drop_diploid = FALSE) {
  out <- facets_tsv(file) |>
    require_cols(c("sample", "gene", "chrom", "gene_start", "tcn", "lcn", "cn_state", "filter"),
      what = basename(file)) |>
    mutate(chrom = normalize_chrom(chrom))
  if (drop_diploid) out <- filter(out, cn_state != "DIPLOID")
  bind_cols(split_pair_name(out$sample), out)
}

#' Read a FACETS arm-level table
#'
#' @param file `*.arm_level.txt`.
#' @return Tibble with `Sample`, `NID`, `arm`, `chrom`, `tcn`, `lcn`,
#'   `frac_of_arm`, `cn_state`.
read_facets_arm <- function(file) {
  out <- facets_tsv(file) |>
    require_cols(c("sample", "arm", "tcn", "lcn", "frac_of_arm", "cn_state"), what = basename(file)) |>
    mutate(
      arm = str_replace(arm, "^23", "X"),
      chrom = normalize_chrom(str_remove(arm, "[pq]$"))
    )
  bind_cols(split_pair_name(out$sample), out) |>
    select(Sample, NID, arm, chrom, everything(), -sample)
}

#' Read a FACETS per-pair QC table
#'
#' @param file `*.facets_qc.txt`.
#' @return One-row tibble: Sample, NID, facets_qc, purity, ploidy, dipLogR,
#'   wgd, fga, n_segs (when present) and the remaining QC columns.
read_facets_qc <- function(file) {
  out <- facets_tsv(file) |>
    require_cols(c("tumor_sample_id", "facets_qc", "purity", "ploidy", "dipLogR"), what = basename(file))
  bind_cols(split_pair_name(out$tumor_sample_id), out) |>
    mutate(facets_qc = as.logical(facets_qc)) |>
    select(Sample, NID, facets_qc, purity, ploidy, dipLogR, any_of(c("wgd", "fga", "n_segs")),
      everything(), -tumor_sample_id, -any_of(c("path", "purity_run_prefix", "hisens_run_prefix")))
}

#' Read a FACETS cncf segment table
#'
#' @param file `*_hisens.cncf.txt` or `*_purity.cncf.txt`.
#' @return Tibble with `Sample`, `chrom` normalized, `loc.start`, `loc.end`,
#'   `tcn`, `lcn`, `cf`, `cnlr.median`, `mafR`, `num.mark`.
read_facets_cncf <- function(file) {
  out <- facets_tsv(file) |>
    require_cols(c("ID", "chrom", "loc.start", "loc.end", "tcn", "lcn"), what = basename(file)) |>
    mutate(chrom = normalize_chrom(chrom))
  ids <- split_pair_name(str_remove(out$ID, "_(hisens|purity)$"))
  bind_cols(ids, out)
}

#' Read a FACETS seg file
#'
#' @param file `*_hisens.seg`.
#' @return Tibble: Sample, NID, chrom, loc.start, loc.end, num.mark, seg.mean.
read_facets_seg <- function(file) {
  out <- facets_tsv(file) |>
    require_cols(c("ID", "chrom", "loc.start", "loc.end", "seg.mean"), what = basename(file)) |>
    mutate(chrom = normalize_chrom(chrom))
  bind_cols(split_pair_name(out$ID), out) |> select(-ID)
}

#' Read the FACETS gene, arm and QC tables for a cohort
#'
#' @param cohort Tibble from `build_cohort()`.
#' @param cache_dir Passed to `cached_read()`.
#' @param drop_diploid Passed to `read_facets_gene()`.
#' @return List with `gene`, `arm`, `qc` tibbles. Sample ids come from the
#'   cohort, not the FACETS `sample` column.
read_cohort_facets <- function(cohort, cache_dir = NULL, drop_diploid = TRUE) {
  files <- cohort_pair_files(cohort, "snv_dir")
  as_paths <- function(col) set_names(files[[col]], files$Sample)
  strip <- function(x) select(x, -any_of(c("Sample_file", "NID_file")))
  gene <- read_cohort_files(cohort, as_paths("facets_gene"), read_facets_gene,
    cache_dir = cache_dir, what = "FACETS gene", drop_diploid = drop_diploid)
  arm <- read_cohort_files(cohort, as_paths("facets_arm"), read_facets_arm,
    cache_dir = cache_dir, what = "FACETS arm")
  qc <- read_cohort_files(cohort, as_paths("facets_qc"), read_facets_qc,
    cache_dir = cache_dir, what = "FACETS qc")
  list(gene = strip(gene), arm = strip(arm), qc = strip(qc))
}
