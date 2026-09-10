# 24_burden.R
#
# Per-sample summary: mutation and SV burden, FACETS purity/ploidy/WGD/FGA,
# the FACETS QC gate and the Tempo QC status.

#' Per-sample summary table
#'
#' @param cohort Tibble from `build_cohort()`.
#' @param snv Tibble from `snv_events()`; `NULL` skips the SNV columns.
#' @param sv Tibble from `sv_events()`; `NULL` skips SV columns.
#' @param cnv_gene Tibble from `cnv_events()`; `NULL` skips CNV columns.
#' @param facets_qc Gated QC table from `facets_qc_gate()`.
#' @param sample_data Tibble from `read_cohort_qc()$sample_data` (TMB, WGD).
#' @param qc_status Tibble from `read_cohort_qc()$qc_status`.
#' @param blank_failed_qc Set purity, ploidy and fga to `NA` when the FACETS
#'   `facets_qc` flag is `FALSE`.
#' @param autosomes_only Count CNV genes on chromosomes 1 to 22 only
#'   (default; see `filter_autosomes()`).
#' @param genome Genome name.
#' @return One row per cohort sample.
sample_summary <- function(cohort, snv = NULL, sv = NULL, cnv_gene = NULL, facets_qc = NULL,
                           sample_data = NULL, qc_status = NULL, blank_failed_qc = TRUE,
                           autosomes_only = TRUE, genome = "hg19") {
  out <- cohort |> select(Sample, NID, ProjNo, has_snv, has_sv, has_facets)

  if (!is.null(snv)) {
    s <- snv |> group_by(Sample) |>
      summarize(n_snv = n(), n_snv_genes = n_distinct(Gene), n_snv_lof = sum(is_lof), .groups = "drop")
    out <- out |> left_join(s, by = "Sample") |>
      mutate(across(c(n_snv, n_snv_genes, n_snv_lof), \(x) if_else(has_snv, coalesce(x, 0L), NA_integer_)))
  }
  if (!is.null(sample_data) && nrow(sample_data) > 0) {
    sd <- sample_data |> select(Sample, any_of(c("TMB", "Number_of_Mutations", "WGD_status", "MSIscore"))) |>
      distinct(Sample, .keep_all = TRUE)
    out <- left_join(out, sd, by = "Sample")
  }
  if (!is.null(sv)) {
    classes <- c("TRA", "DEL", "DUP", "INV", "INS")
    v <- sv |> distinct(Sample, UUID, sv_class, is_fusion) |>
      group_by(Sample) |>
      summarize(
        n_sv = n(),
        n_sv_TRA = sum(sv_class == "TRA"), n_sv_DEL = sum(sv_class == "DEL"),
        n_sv_DUP = sum(sv_class == "DUP"), n_sv_INV = sum(sv_class == "INV"),
        n_sv_INS = sum(sv_class == "INS"), n_fusion = sum(is_fusion),
        .groups = "drop"
      )
    out <- out |> left_join(v, by = "Sample") |>
      mutate(across(c(n_sv, starts_with("n_sv_"), n_fusion), \(x) if_else(has_sv, coalesce(x, 0L), NA_integer_)))
  }
  if (!is.null(facets_qc) && nrow(facets_qc) > 0) {
    fq <- facets_qc |>
      select(Sample, facets_qc, purity, ploidy, dipLogR, any_of(c("wgd", "fga", "n_segs", "EXCLUDE", "exclude_reason"))) |>
      distinct(Sample, .keep_all = TRUE)
    if (blank_failed_qc) {
      fq <- fq |> mutate(across(any_of(c("purity", "ploidy", "fga")),
        \(x) if_else(coalesce(facets_qc, FALSE), x, NA_real_)))
    }
    out <- left_join(out, fq, by = "Sample")
  }
  if (!is.null(cnv_gene) && nrow(cnv_gene) > 0) {
    if (autosomes_only) cnv_gene <- filter_autosomes(cnv_gene, "chrom", genome)
    cn <- cnv_gene |> distinct(Sample, Gene, Dir) |>
      group_by(Sample) |>
      summarize(n_cnv_genes = n_distinct(Gene),
        n_cnv_gain = n_distinct(Gene[Dir %in% "gain"]),
        n_cnv_loss = n_distinct(Gene[Dir %in% "loss"]), .groups = "drop")
    out <- left_join(out, cn, by = "Sample")
  }
  if (!is.null(qc_status) && nrow(qc_status) > 0) {
    qs <- qc_status
    if ("is_tumor" %in% names(qs)) qs <- filter(qs, is_tumor)
    qs <- qs |> filter(Sample %in% cohort$Sample) |>
      distinct(Sample, .keep_all = TRUE) |>
      select(Sample, QC_Status = Status, QC_Reason = Reason)
    out <- left_join(out, qs, by = "Sample")
  }
  arrange(out, ProjNo, Sample)
}
