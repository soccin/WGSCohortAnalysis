# 02_tempo_paths.R
#
# Locate the files of one Tempo tumor/normal pair from any one of its paths.
# Replaces the five discovery idioms in the older scripts
# (dirname(dirname(maf))/facets, nested dir_ls, zip and tar walks).
#
# Per-pair layout (both vintages seen so far):
#
#   <run>/somatic/<PAIR>/
#     combined_mutations/<PAIR>.somatic.final.maf     (+ .somatic.unfiltered.maf)
#     combined_svs/<PAIR>.final.bedpe
#     facets/<PAIR>/<PAIR>.facets_qc.txt
#     facets/<PAIR>/<PAIR>_OUT.txt
#     facets/<PAIR>/facets<ver>/<PAIR>.gene_level.txt
#     facets/<PAIR>/facets<ver>/<PAIR>.arm_level.txt
#     facets/<PAIR>/facets<ver>/<PAIR>_hisens.cncf.txt   (+ _purity)
#     facets/<PAIR>/facets<ver>/<PAIR>_hisens.seg
#     facets/<PAIR>/facets<ver>/<PAIR>.qc.txt
#     meta_data/<PAIR>.sample_data.txt                 (newer runs)
#     multiqc/<PAIR>.QC_Status.txt                     (older runs)
#   <run>/cohort_level/<cohort>/sample_data.txt, alignment_qc.txt

#' Resolve the per-pair directory from a path inside it
#'
#' @param path A MAF, bedpe, or the pair directory itself.
#' @return The `somatic/<PAIR>` directory.
tempo_pair_dir <- function(path) {
  path <- as.character(path)
  if (fs::is_dir(path)) return(path)
  as.character(fs::path_dir(fs::path_dir(path)))
}

#' Expected files of a Tempo pair
#'
#' Every element is a path or `NA` when the file is absent. `facets_*`
#' searches are recursive under `facets/` so the versioned subdirectory name
#' does not matter. When two fits exist, `hisens` is preferred for cncf/seg.
#'
#' @param path A MAF, bedpe, or pair directory.
#' @return A one-row tibble: pair_dir, pair, Sample, NID, maf, maf_unfiltered,
#'   bedpe, facets_gene, facets_arm, facets_qc, facets_out, facets_cncf,
#'   facets_seg, facets_fit_qc, sample_data, qc_status.
tempo_pair_files <- function(path) {
  dir <- tempo_pair_dir(path)
  pair <- fs::path_file(dir)
  first_or_na <- function(sub, regexp, recurse = FALSE) {
    d <- fs::path(dir, sub)
    if (!fs::dir_exists(d)) return(NA_character_)
    hits <- fs::dir_ls(d, regexp = regexp, recurse = recurse, type = "file")
    if (length(hits) == 0) NA_character_ else as.character(sort(hits)[[1]])
  }
  ids <- split_pair_name(pair)
  tibble(
    pair_dir = dir,
    pair = pair,
    Sample = ids$Sample,
    NID = ids$NID,
    maf = first_or_na("combined_mutations", "\\.somatic\\.final\\.maf$"),
    maf_unfiltered = first_or_na("combined_mutations", "\\.somatic\\.unfiltered\\.maf$"),
    bedpe = first_or_na("combined_svs", "\\.final\\.bedpe$"),
    facets_gene = first_or_na("facets", "\\.gene_level\\.txt$", recurse = TRUE),
    facets_arm = first_or_na("facets", "\\.arm_level\\.txt$", recurse = TRUE),
    facets_qc = first_or_na("facets", "\\.facets_qc\\.txt$", recurse = TRUE),
    facets_out = first_or_na("facets", "_OUT\\.txt$", recurse = TRUE),
    facets_cncf = first_or_na("facets", "_hisens\\.cncf\\.txt$", recurse = TRUE),
    facets_seg = first_or_na("facets", "_hisens\\.seg$", recurse = TRUE),
    facets_fit_qc = first_or_na("facets", "[^_]\\.qc\\.txt$", recurse = TRUE),
    sample_data = first_or_na("meta_data", "\\.sample_data\\.txt$"),
    qc_status = first_or_na("multiqc", "\\.QC_Status\\.txt$")
  )
}

#' Cohort-level files of the Tempo run that contains a pair
#'
#' @param path A MAF, bedpe, or pair directory.
#' @return One-row tibble: run_dir, sample_data, alignment_qc (paths or `NA`).
tempo_cohort_files <- function(path) {
  run_dir <- fs::path_dir(fs::path_dir(tempo_pair_dir(path)))
  find1 <- function(regexp) {
    d <- fs::path(run_dir, "cohort_level")
    if (!fs::dir_exists(d)) return(NA_character_)
    hits <- fs::dir_ls(d, regexp = regexp, recurse = TRUE, type = "file")
    if (length(hits) == 0) NA_character_ else as.character(sort(hits)[[1]])
  }
  tibble(
    run_dir = as.character(run_dir),
    sample_data = find1("/sample_data\\.txt$"),
    alignment_qc = find1("/alignment_qc\\.txt$")
  )
}

#' Pair files for every row of a cohort
#'
#' @param cohort Tibble from `build_cohort()`.
#' @param side Which directory to resolve from: `snv_dir` (default, holds
#'   MAF and FACETS) or `sv_dir`.
#' @return `tempo_pair_files()` rows bound with the cohort `Sample`.
cohort_pair_files <- function(cohort, side = c("snv_dir", "sv_dir")) {
  side <- match.arg(side)
  cohort |>
    filter(!is.na(.data[[side]])) |>
    select(Sample, dir = all_of(side)) |>
    mutate(files = map(dir, \(d) tempo_pair_files(d) |> select(-Sample, -NID))) |>
    select(Sample, files) |>
    tidyr::unnest(files)
}
