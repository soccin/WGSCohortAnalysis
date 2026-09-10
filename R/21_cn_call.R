# 21_cn_call.R
#
# Copy-number direction calls and the FACETS sample QC gate. The five
# incompatible amp/del rules in the older scripts are captured as presets so
# a report states which one it used.

#' Copy-number direction from FACETS columns
#'
#' Presets:
#' * `facets` (default): from `cn_state` with the regexes used in the AML
#'   OGM comparison. loss = HOMDEL, HETLOSS, LOSS..., DOUBLE LOSS...;
#'   gain = AMP..., GAIN...; cnloh = CNLOH. Loss wins, then gain.
#' * `strict`: the 260719 report rule. gain = GAIN only, loss = HETLOSS or
#'   HOMDEL only; AMP, CNLOH and the WGD-relative states are `NA`.
#' * `tcn`: loss = tcn < 2, gain = tcn > 2, cnloh = tcn == 2 & lcn == 0.
#' * `focal`: the `tcn` rule restricted to PASS/RESCUE calls on segments
#'   with fewer than `max_genes_on_seg` genes; everything else `NA`.
#'
#' @param cn_state FACETS `cn_state`.
#' @param tcn,lcn Total and lesser copy number.
#' @param filter FACETS `filter` column (needed by `focal`).
#' @param genes_on_seg Needed by `focal`.
#' @param preset One of the presets above.
#' @param max_genes_on_seg Focal threshold.
#' @return Character vector: `"gain"`, `"loss"`, `"cnloh"` or `NA`.
cn_call <- function(cn_state, tcn = NULL, lcn = NULL, filter = NULL, genes_on_seg = NULL,
                    preset = c("facets", "strict", "tcn", "focal"), max_genes_on_seg = 10) {
  preset <- match.arg(preset)
  if (preset == "facets") {
    return(case_when(
      str_detect(cn_state, "HOMDEL|HETLOSS|^LOSS|DOUBLE LOSS") ~ "loss",
      str_detect(cn_state, "^AMP|^GAIN|GAIN$") ~ "gain",
      str_detect(cn_state, "CNLOH") ~ "cnloh",
      TRUE ~ NA_character_
    ))
  }
  if (preset == "strict") {
    return(case_when(
      cn_state == "GAIN" ~ "gain",
      cn_state %in% c("HETLOSS", "HOMDEL") ~ "loss",
      TRUE ~ NA_character_
    ))
  }
  if (is.null(tcn)) wca_abort("cn_call(): preset '{preset}' needs tcn")
  tcn <- suppressWarnings(as.numeric(tcn))
  lcn <- if (is.null(lcn)) rep(NA_real_, length(tcn)) else suppressWarnings(as.numeric(lcn))
  dir <- case_when(
    tcn < 2 ~ "loss",
    tcn > 2 ~ "gain",
    tcn == 2 & !is.na(lcn) & lcn == 0 ~ "cnloh",
    TRUE ~ NA_character_
  )
  if (preset == "focal") {
    if (is.null(filter) || is.null(genes_on_seg)) {
      wca_abort("cn_call(): preset 'focal' needs filter and genes_on_seg")
    }
    ok <- filter %in% c("PASS", "RESCUE") &
      coalesce(suppressWarnings(as.numeric(genes_on_seg)) < max_genes_on_seg, FALSE)
    dir <- if_else(ok, dir, NA_character_)
  }
  dir
}

#' FACETS sample QC gate
#'
#' Adds `EXCLUDE` and `exclude_reason` to a FACETS QC table. A sample is
#' excluded when `abs(dipLogR) > max_abs_diplogr` (unreliable ploidy) and,
#' optionally, when the FACETS `facets_qc` flag is `FALSE`.
#'
#' @param facets_qc Tibble from `read_facets_qc()` / `read_cohort_facets()`.
#' @param max_abs_diplogr dipLogR cutoff; `NULL` or `Inf` disables.
#' @param require_facets_qc Also exclude `facets_qc == FALSE`.
#' @return The table with `EXCLUDE` (logical) and `exclude_reason`.
facets_qc_gate <- function(facets_qc, max_abs_diplogr = 1.5, require_facets_qc = FALSE) {
  cut <- max_abs_diplogr %||% Inf
  facets_qc |>
    mutate(
      .dip = coalesce(abs(dipLogR) > cut, FALSE),
      .qc = require_facets_qc & coalesce(!facets_qc, TRUE),
      EXCLUDE = .dip | .qc,
      exclude_reason = case_when(
        .dip & .qc ~ str_glue("|dipLogR| > {cut}; facets_qc FALSE"),
        .dip ~ str_glue("|dipLogR| > {cut}"),
        .qc ~ "facets_qc FALSE",
        TRUE ~ NA_character_
      ) |> as.character()
    ) |>
    select(-.dip, -.qc)
}

#' Samples passing the FACETS gate
#'
#' @param facets_qc Table with an `EXCLUDE` column (see `facets_qc_gate()`).
#' @return Character vector of `Sample`.
cnv_samples <- function(facets_qc) {
  require_cols(facets_qc, c("Sample", "EXCLUDE"), "facets_qc")
  facets_qc |> filter(!EXCLUDE) |> pull(Sample) |> unique()
}
