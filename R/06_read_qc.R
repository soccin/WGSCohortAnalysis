# 06_read_qc.R
#
# Tempo QC readers: the per-pair QC_Status.txt (older runs), the per-pair
# meta_data sample_data.txt (newer runs; TMB, WGD, MSI) and the cohort-level
# alignment_qc.txt.

#' Read a Tempo QC_Status file
#'
#' The file has an unnamed first column holding the sample id and one row per
#' tumor and normal.
#'
#' @param file `multiqc/<PAIR>.QC_Status.txt`.
#' @return Tibble: Sample, Status, Reason.
read_qc_status <- function(file) {
  header <- str_split_1(readr::read_lines(file, n_max = 1), "\t")
  header[header == ""] <- "Sample"
  read_tsv(file, skip = 1, col_names = header, col_types = cols(.default = "c"), progress = FALSE) |>
    require_cols(c("Sample", "Status", "Reason"), what = basename(file)) |>
    select(Sample, Status, Reason)
}

#' Read a Tempo sample_data table
#'
#' Works for the per-pair `meta_data/<PAIR>.sample_data.txt` and the
#' cohort-level `sample_data.txt`. Signature (`SBS*`) and HLA columns are
#' dropped unless `keep_all = TRUE`.
#'
#' @param file Path.
#' @param keep_all Keep every column.
#' @return Tibble: Sample, NID, purity, ploidy, WGD_status, MSIscore,
#'   Number_of_Mutations, TMB (when present).
read_sample_data <- function(file, keep_all = FALSE) {
  out <- read_tsv(file, col_types = cols(.default = "c"), progress = FALSE) |>
    require_cols("sample", what = basename(file)) |>
    mutate(across(everything(), str_trim)) |>
    type_convert_silent()
  if ("WGD_status" %in% names(out)) {
    out <- mutate(out, WGD_status = as.logical(toupper(as.character(WGD_status))))
  }
  if (!keep_all) out <- select(out, -matches("^SBS"), -matches("^HLA"))
  bind_cols(split_pair_name(out$sample), out) |>
    select(Sample, NID, any_of(c("purity", "ploidy", "WGD_status", "MSIscore",
      "MSI_Total_Sites", "MSI_Somatic_Sites", "Number_of_Mutations", "TMB")),
      everything(), -sample)
}

#' Read a Tempo cohort-level alignment_qc table
#'
#' @param file `cohort_level/<cohort>/alignment_qc.txt`.
#' @return Typed tibble, one row per sample (tumors and normals).
read_alignment_qc <- function(file) {
  read_tsv(file, col_types = cols(.default = "c"), progress = FALSE) |>
    require_cols(c("Sample", "MedianCoverage"), what = basename(file)) |>
    type_convert_silent()
}

#' Read QC_Status and sample_data for every cohort pair
#'
#' @param cohort Tibble from `build_cohort()`.
#' @param cache_dir Passed to `cached_read()`.
#' @return List with `qc_status` (tumor and normal rows; `is_tumor` marks
#'   the tumor row, `id` is the id in the file) and `sample_data`. Either can be an empty tibble when the run vintage does
#'   not emit the file.
read_cohort_qc <- function(cohort, cache_dir = NULL) {
  files <- cohort_pair_files(cohort, "snv_dir")
  as_paths <- function(col) set_names(files[[col]], files$Sample)
  qc_status <- read_cohort_files(cohort, as_paths("qc_status"), read_qc_status,
    cache_dir = cache_dir, what = "QC_Status")
  if (nrow(qc_status) > 0) {
    qc_status <- qc_status |>
      mutate(is_tumor = Sample_file == Sample) |>
      select(Sample, NID, id = Sample_file, is_tumor, Status, Reason)
  }
  sample_data <- read_cohort_files(cohort, as_paths("sample_data"), read_sample_data,
    cache_dir = cache_dir, what = "sample_data")
  sample_data <- select(sample_data, -any_of(c("Sample_file", "NID_file")))
  list(qc_status = qc_status, sample_data = sample_data)
}
