# 04_read_sv.R
#
# Tempo somatic SV (bedpe) reader. Consolidates the six older variants: the
# INFO_A key=value block and the caller-specific FORMAT fields (t_ / n_ for
# tumor / normal) are pivoted into columns, VAFs are derived per caller, and
# the SV class comes from the ID prefix so BND rows read as TRA.

#' Priority column order for SV tables
#'
#' @return Character vector read from `data/sv_report_cols.txt`.
sv_report_cols <- function() {
  read_list_file(wca_file("data", "sv_report_cols.txt"))
}

#' SV class from a Tempo ID
#'
#' IDs look like `TEMPO_DEL_1_804712_1_804990_+-`; the second token is one of
#' DEL, DUP, INV, INS, TRA. The `TYPE` column uses BND for translocations, so
#' this is the column to summarize on.
#'
#' @param id Character vector of IDs.
#' @return Character vector of classes.
sv_class <- function(id) {
  str_split_i(id, "_", 2)
}

#' Read one Tempo somatic bedpe
#'
#' @param file Path to `*.final.bedpe`.
#' @param detail `"full"` keeps every parsed column; `"standard"` keeps the
#'   report columns, VAFs and read support; `"minimal"` keeps only the
#'   priority columns plus support counts.
#' @param priority_cols Columns placed first, in this order.
#' @param genome Used to normalize chromosome names.
#' @return Typed tibble with one row per SV, or `NULL` for a file with no
#'   calls. Adds `UUID` (TUMOR_ID_ID), `sv_class`, `GenePair`, `is_inter`
#'   and `span`.
read_tempo_sv <- function(file, detail = c("full", "standard", "minimal"),
                          priority_cols = sv_report_cols(), genome = "hg19") {
  detail <- match.arg(detail)
  head_lines <- readr::read_lines(file, n_max = 500)
  n_skip <- sum(str_starts(head_lines, "##"))
  tt <- read_tsv(file, skip = n_skip, col_types = cols(.default = "c"),
    na = c("", "NA"), progress = FALSE)
  if (nrow(tt) == 0) return(NULL)
  tt <- tt |>
    rename(CHROM_A = `#CHROM_A`) |>
    rename_with(\(x) str_replace(x, "^repName-repClass-repFamily:-", "repeat."), matches("^repName")) |>
    mutate(UUID = str_c(TUMOR_ID, ID, sep = "_")) |>
    select(CHROM_A:END_B, TYPE, UUID, everything())

  info_a <- tt |>
    select(UUID, INFO_A) |>
    tidyr::separate_longer_delim(INFO_A, delim = ";") |>
    filter(str_detect(INFO_A, "=")) |>
    tidyr::separate_wider_delim(INFO_A, delim = "=", names = c("key", "val"),
      too_many = "merge") |>
    distinct(UUID, key, .keep_all = TRUE) |>
    tidyr::pivot_wider(names_from = key, values_from = val)

  pivot_format <- function(col, prefix) {
    fmt <- tt |>
      select(UUID, FORMAT) |>
      tidyr::separate_longer_delim(FORMAT, delim = ":") |>
      mutate(RID = row_number())
    vals <- tt |>
      select(UUID, val = all_of(col)) |>
      tidyr::separate_longer_delim(val, delim = ":") |>
      mutate(RID = row_number())
    full_join(fmt, vals, by = join_by(UUID, RID)) |>
      select(-RID) |>
      mutate(FORMAT = str_c(prefix, FORMAT)) |>
      tidyr::pivot_wider(names_from = FORMAT, values_from = val)
  }
  tumor_fmt <- pivot_format("TUMOR", "t_")
  normal_fmt <- pivot_format("NORMAL", "n_")

  sv <- tt |>
    left_join(info_a, by = "UUID") |>
    left_join(tumor_fmt, by = "UUID") |>
    left_join(normal_fmt, by = "UUID") |>
    select(any_of(priority_cols), everything()) |>
    select(-any_of(c("INFO_A", "INFO_B", "FORMAT", "TUMOR", "NORMAL"))) |>
    mutate(across(c(CHROM_A, CHROM_B), normalize_chrom)) |>
    type_convert_silent() |>
    sv_vaf() |>
    mutate(
      sv_class = sv_class(ID),
      GenePair = str_c(gene1, "::", gene2),
      is_inter = CHROM_A != CHROM_B,
      span = if_else(is_inter, NA_real_, abs(as.numeric(START_B) - as.numeric(START_A)))
    ) |>
    select(GenePair, any_of(priority_cols), sv_class, is_inter, span, everything())

  switch(detail,
    full = sv,
    standard = sv |>
      select(GenePair, any_of(priority_cols), sv_class, is_inter, span, ID,
        matches("VAF$"), matches("^[tn]_.*_(AD|PE|SR|PR|PS|DR|DV|RR|RV|DP)$"),
        matches("^CC_|^DGv|^repeat\\."), any_of(c("Cosmic_Fusion_Counts", "transcript1", "transcript2", "STRAND_A", "STRAND_B", "FILTER"))),
    minimal = sv |>
      select(GenePair, any_of(priority_cols), sv_class, is_inter, span, ID, matches("VAF$"))
  )
}

#' Add per-caller VAF columns to a parsed SV table
#'
#' Delly: span (DV/(DV+DR)) and junction (RV/(RV+RR)) VAFs for tumor and
#' normal. Manta: junction VAF from the `REF,ALT` split-read pair. SvABA:
#' AD/DP. Columns are only computed when the inputs exist, so older bedpes
#' without a caller still parse.
#'
#' @param sv Tibble from `read_tempo_sv()` before VAFs.
#' @return The tibble with `*VAF` columns added.
sv_vaf <- function(sv) {
  has <- function(...) all(c(...) %in% names(sv))
  num <- function(x) suppressWarnings(as.numeric(x))
  if (has("t_delly_DV", "t_delly_DR")) {
    sv <- mutate(sv, t_delly_SpanVAF = num(t_delly_DV) / (num(t_delly_DV) + num(t_delly_DR)))
  }
  if (has("t_delly_RV", "t_delly_RR")) {
    sv <- mutate(sv, t_delly_JuncVAF = num(t_delly_RV) / (num(t_delly_RV) + num(t_delly_RR)))
  }
  if (has("n_delly_DV", "n_delly_DR")) {
    sv <- mutate(sv, n_delly_SpanVAF = num(n_delly_DV) / (num(n_delly_DV) + num(n_delly_DR)))
  }
  if (has("n_delly_RV", "n_delly_RR")) {
    sv <- mutate(sv, n_delly_JuncVAF = num(n_delly_RV) / (num(n_delly_RV) + num(n_delly_RR)))
  }
  if (has("t_manta_SR")) {
    sv <- sv |>
      mutate(
        t_manta_SRR = num(str_split_i(t_manta_SR, ",", 1)),
        t_manta_SRV = num(str_split_i(t_manta_SR, ",", 2)),
        t_manta_JuncVAF = t_manta_SRV / (t_manta_SRV + t_manta_SRR)
      )
  }
  if (has("t_manta_PR")) {
    sv <- sv |>
      mutate(
        t_manta_PRR = num(str_split_i(t_manta_PR, ",", 1)),
        t_manta_PRV = num(str_split_i(t_manta_PR, ",", 2)),
        t_manta_SpanVAF = t_manta_PRV / (t_manta_PRV + t_manta_PRR)
      )
  }
  if (has("t_svaba_AD", "t_svaba_DP")) {
    sv <- mutate(sv, t_svaba_VAF = num(t_svaba_AD) / num(t_svaba_DP))
  }
  if (has("n_svaba_AD", "n_svaba_DP")) {
    sv <- mutate(sv, n_svaba_VAF = num(n_svaba_AD) / num(n_svaba_DP))
  }
  sv
}

#' Cross-caller read support summaries
#'
#' Callers count the same reads, so support is never summed. Adds the
#' MAX and MEDIAN of paired and split read counts across Delly, Manta and
#' SvABA, plus evidence flags.
#'
#' @param sv Tibble from `read_tempo_sv()`.
#' @param min_split,min_paired Thresholds for the `meets_*` flags.
#' @return The tibble with `max_paired_reads`, `max_split_reads`,
#'   `median_paired_reads`, `median_split_reads`, `total_max_reads`,
#'   `has_paired_evidence`, `has_split_evidence`, `multi_evidence_type`,
#'   `meets_split_threshold`, `meets_paired_threshold`.
sv_support <- function(sv, min_split = 3, min_paired = 5) {
  num_col <- function(name) {
    if (name %in% names(sv)) suppressWarnings(as.numeric(sv[[name]])) else rep(NA_real_, nrow(sv))
  }
  paired <- cbind(num_col("t_delly_DV"), num_col("t_manta_PRV"), num_col("t_svaba_DR"))
  split <- cbind(num_col("t_delly_RV"), num_col("t_manta_SRV"), num_col("t_svaba_SR"))
  row_stat <- function(m, f) {
    out <- apply(m, 1, \(x) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE))
    as.numeric(out)
  }
  sv |>
    mutate(
      max_paired_reads = row_stat(paired, max),
      max_split_reads = row_stat(split, max),
      median_paired_reads = row_stat(paired, median),
      median_split_reads = row_stat(split, median),
      total_max_reads = coalesce(max_paired_reads, 0) + coalesce(max_split_reads, 0),
      has_paired_evidence = coalesce(max_paired_reads > 0, FALSE),
      has_split_evidence = coalesce(max_split_reads > 0, FALSE),
      multi_evidence_type = has_paired_evidence & has_split_evidence,
      meets_split_threshold = coalesce(max_split_reads >= min_split, FALSE),
      meets_paired_threshold = coalesce(max_paired_reads >= min_paired, FALSE)
    )
}

#' Read the bedpe of every cohort pair
#'
#' @inheritParams read_cohort_maf
#' @param ... Passed to `read_tempo_sv()`.
#' @return Bound tibble with `Sample`, `NID` first. Pairs whose bedpe has no
#'   calls contribute no rows.
read_cohort_sv <- function(cohort, cache_dir = NULL, ...) {
  read_cohort_files(cohort, file_col = "sv_file", reader = read_tempo_sv,
    cache_dir = cache_dir, what = "SV", ...)
}
