# 23_recurrence.R
#
# Per-gene, per-band, per-arm and per-pair recurrence. Every table carries
# its denominator as an attribute so a report can state it, and PCT is
# always n / denominator with n counting samples, not events.

#' Per-gene recurrence
#'
#' @param tbl Table with `Sample` and a gene column.
#' @param denominator Number of samples the frequency is relative to.
#' @param event Label stored in the `Event` column.
#' @param gene_col Gene column name.
#' @return Tibble: Gene, n, PCT, Event, sorted by n. Attribute
#'   `denominator` holds the denominator.
recurrence_table <- function(tbl, denominator, event, gene_col = "Gene") {
  if (missing(denominator) || is.null(denominator)) wca_abort("recurrence_table(): denominator is required")
  out <- tbl |>
    distinct(Sample, Gene = .data[[gene_col]]) |>
    count(Gene) |>
    mutate(PCT = n / denominator, Event = event) |>
    arrange(desc(n), Gene)
  attr(out, "denominator") <- denominator
  out
}

#' Wide per-gene summary across event types
#'
#' One row per gene with `<Event>_n` and `<Event>_PCT` for every table in
#' `tables`, `Nt` = sum of the counts and `PCT = Nt / n_samples`. As in the
#' original reports `Nt` sums event types, so it can exceed the sample count
#' when a gene carries several event types in one sample.
#'
#' @param tables Named list of `recurrence_table()` results.
#' @param n_samples Cohort size for the overall PCT.
#' @return Tibble: Gene, Nt, PCT, then `*_n` and `*_PCT` columns.
event_summary <- function(tables, n_samples) {
  tables <- purrr::compact(tables)
  if (length(tables) == 0) return(tibble(Gene = character(), Nt = numeric(), PCT = numeric()))
  long <- bind_rows(tables) |>
    tidyr::pivot_longer(c(n, PCT), names_to = "K", values_to = "V") |>
    tidyr::unite("EK", Event, K, sep = "_")
  wide <- long |>
    tidyr::pivot_wider(names_from = EK, values_from = V, values_fill = 0)
  n_cols <- sort(str_subset(names(wide), "_n$"))
  pct_cols <- sort(str_subset(names(wide), "_PCT$"))
  wide |>
    mutate(Nt = rowSums(across(all_of(n_cols)), na.rm = TRUE), PCT = Nt / n_samples) |>
    arrange(desc(Nt), Gene) |>
    select(Gene, Nt, PCT, all_of(n_cols), all_of(pct_cols))
}

#' Filter a summary to a gene set
#'
#' @param summary Tibble from `event_summary()`.
#' @param genes Exact gene symbols (not a regex; `DAD1` should not match `ADAD1`).
#' @return Rows for `genes`, alphabetically.
gene_set_table <- function(summary, genes) {
  summary |> filter(Gene %in% genes) |> arrange(Gene)
}

#' Cytoband gain/loss recurrence
#'
#' Each gene is mapped to its band by `gene_start`; a sample counts once per
#' band and direction however many genes it hits. Ranked by the combined
#' count, ties in genome order.
#'
#' @param cnv_gene Tibble from `cnv_events()` (needs `Dir`).
#' @param n_cnv_samples Denominator: FACETS-passing samples.
#' @param genome Genome name.
#' @param autosomes_only Drop X and Y rows before counting (default; see
#'   `filter_autosomes()`).
#' @return Tibble: Band, chrom, N_total, N_gain, PCT_gain, N_loss, PCT_loss
#'   (and N_cnloh, PCT_cnloh when present).
cytoband_recurrence <- function(cnv_gene, n_cnv_samples, genome = "hg19", autosomes_only = TRUE) {
  if (autosomes_only) cnv_gene <- filter_autosomes(cnv_gene, "chrom", genome)
  bands <- cnv_gene |>
    filter(!is.na(Dir)) |>
    annotate_cytoband("chrom", "gene_start", genome, name = "Band", keep_bounds = TRUE) |>
    filter(!is.na(Band)) |>
    distinct(Band, chrom, Band_start, Sample, Dir) |>
    count(Band, chrom, Band_start, Dir) |>
    tidyr::pivot_wider(names_from = Dir, values_from = n, names_prefix = "N_", values_fill = 0)
  for (col in c("N_gain", "N_loss")) if (!col %in% names(bands)) bands[[col]] <- 0L
  out <- bands |>
    mutate(
      N_total = rowSums(across(starts_with("N_"))),
      PCT_gain = N_gain / n_cnv_samples,
      PCT_loss = N_loss / n_cnv_samples
    )
  if ("N_cnloh" %in% names(out)) out <- mutate(out, PCT_cnloh = N_cnloh / n_cnv_samples)
  out <- out |>
    arrange(desc(N_total), chrom_factor(chrom, genome), Band_start) |>
    select(Band, chrom, N_total, N_gain, PCT_gain, N_loss, PCT_loss, any_of(c("N_cnloh", "PCT_cnloh")))
  attr(out, "denominator") <- n_cnv_samples
  out
}

#' Arm-level recurrence
#'
#' An arm counts as altered in a sample when the FACETS arm call covers at
#' least `min_frac` of the arm.
#'
#' @param cnv_arm Tibble from `cnv_arm_events()`.
#' @param n_cnv_samples Denominator.
#' @param min_frac Minimum `frac_of_arm`.
#' @param genome Genome name.
#' @param autosomes_only Restrict both the arm list and the calls to
#'   chromosomes 1 to 22 (default; see `filter_autosomes()`).
#' @return Tibble: arm, chrom, N_gain, PCT_gain, N_loss, PCT_loss
#'   (and cnloh), in genome order.
arm_recurrence <- function(cnv_arm, n_cnv_samples, min_frac = 0.5, genome = "hg19",
                           autosomes_only = TRUE) {
  arms <- chrom_arms(genome) |> select(arm, chrom)
  if (autosomes_only) {
    arms <- filter_autosomes(arms, "chrom", genome)
    cnv_arm <- filter_autosomes(cnv_arm, "chrom", genome)
  }
  calls <- cnv_arm |>
    filter(!is.na(Dir), coalesce(frac_of_arm >= min_frac, FALSE)) |>
    distinct(Sample, arm, Dir) |>
    count(arm, Dir) |>
    tidyr::pivot_wider(names_from = Dir, values_from = n, names_prefix = "N_", values_fill = 0)
  out <- arms |> left_join(calls, by = "arm")
  for (col in c("N_gain", "N_loss")) if (!col %in% names(out)) out[[col]] <- 0L
  out <- out |>
    mutate(across(starts_with("N_"), \(x) coalesce(x, 0L)),
      PCT_gain = N_gain / n_cnv_samples, PCT_loss = N_loss / n_cnv_samples)
  if ("N_cnloh" %in% names(out)) out <- mutate(out, PCT_cnloh = N_cnloh / n_cnv_samples)
  out <- out |>
    arrange(chrom_factor(chrom, genome), arm) |>
    select(arm, chrom, N_gain, PCT_gain, N_loss, PCT_loss, any_of(c("N_cnloh", "PCT_cnloh")))
  attr(out, "denominator") <- n_cnv_samples
  out
}

#' Recurrent gene pairs across samples
#'
#' @param sv Tibble from `sv_events()` (needs `PairKey`).
#' @param n_samples Denominator.
#' @param min_n Keep pairs seen in at least this many samples.
#' @return Tibble: PairKey, gene1, gene2, n, PCT, Classes, Samples.
sv_pair_recurrence <- function(sv, n_samples, min_n = 1) {
  out <- sv |>
    distinct(PairKey, Sample, sv_class) |>
    group_by(PairKey) |>
    summarize(
      n = n_distinct(Sample),
      Classes = str_c(sort(unique(sv_class)), collapse = ";"),
      Samples = str_c(sort(unique(Sample)), collapse = ";"),
      .groups = "drop"
    ) |>
    mutate(PCT = n / n_samples, gene1 = str_split_i(PairKey, "::", 1), gene2 = str_split_i(PairKey, "::", 2)) |>
    filter(n >= min_n) |>
    arrange(desc(n), PairKey) |>
    select(PairKey, gene1, gene2, n, PCT, Classes, Samples)
  attr(out, "denominator") <- n_samples
  out
}

#' SV partner table per gene
#'
#' `N` counts distinct Gene-Sample-Partner triples, not samples, so a gene
#' seen with several partners in one sample can outrank a gene recurrent
#' across samples (kept for continuity with the 260719 report).
#'
#' @param svg Tibble from `sv_gene_events()`.
#' @return Tibble: Gene, N, GenePartners, Partners.
sv_partner_table <- function(svg) {
  svg |>
    distinct(Gene, Sample, gene1, gene2) |>
    tidyr::pivot_longer(c(gene1, gene2), names_to = "Side", values_to = "GeneA") |>
    filter(Gene != GeneA) |>
    distinct(Gene, Sample, GeneA) |>
    group_by(Gene) |>
    summarize(
      N = n(),
      GenePartners = str_c(unique(sort(GeneA)), collapse = ";"),
      Partners = str_c(str_c(GeneA, "::", Sample), collapse = ";"),
      .groups = "drop"
    ) |>
    arrange(desc(N), Gene)
}

#' Fusion table
#'
#' @param sv Tibble from `sv_events()`.
#' @param n_samples Denominator.
#' @param classes Fusion classes to keep.
#' @return Tibble: GenePair, fusion_class, n, PCT, Samples, Classes.
fusion_table <- function(sv, n_samples, classes = c("in-frame", "out-of-frame", "protein-fusion", "transcript")) {
  out <- sv |>
    filter(fusion_class %in% classes) |>
    distinct(GenePair, fusion_class, Sample, sv_class) |>
    group_by(GenePair, fusion_class) |>
    summarize(
      n = n_distinct(Sample),
      Classes = str_c(sort(unique(sv_class)), collapse = ";"),
      Samples = str_c(sort(unique(Sample)), collapse = ";"),
      .groups = "drop"
    ) |>
    mutate(PCT = n / n_samples) |>
    arrange(desc(n), GenePair) |>
    select(GenePair, fusion_class, n, PCT, Classes, Samples)
  attr(out, "denominator") <- n_samples
  out
}

#' Recurrence split by sample group
#'
#' Counts samples per key (a gene, a `PairKey`, a `GenePair`) within each
#' named group of samples, so a question such as "hit in every cell line
#' and in at least one patient" is a filter on the result. A sample counts
#' once per key however many events it carries; samples outside every
#' group are ignored, and `n` counts grouped samples only.
#'
#' @param tbl Table with `Sample` and `key_col`.
#' @param groups Named list, group name to character vector of samples.
#'   Group names become column suffixes, so keep them simple (`CellLine`,
#'   `Patient`).
#' @param key_col Column to count by.
#' @param detail_cols Columns whose distinct values are collapsed with `;`
#'   into a column of the same name (e.g. `sv_class`, `fusion_class`).
#' @return Tibble: `<key_col>`, n, then `n_<group>`, `PCT_<group>` and
#'   `Samples_<group>` per group, then `detail_cols`, sorted by n and by the
#'   group counts in the order given. Attribute `denominator` is the named
#'   vector of group sizes.
group_recurrence <- function(tbl, groups, key_col = "Gene", detail_cols = character()) {
  if (is.null(names(groups)) || any(names(groups) == "")) wca_abort("group_recurrence(): groups must be named")
  if (any(str_detect(names(groups), "[^A-Za-z0-9]"))) wca_abort("group_recurrence(): group names must be alphanumeric")
  require_cols(tbl, c("Sample", key_col, detail_cols), "tbl")
  sizes <- purrr::map_int(groups, \(x) length(unique(x)))
  membership <- purrr::imap(groups, \(x, g) tibble(Sample = unique(x), Group = g)) |> bind_rows()
  dup <- membership |> count(Sample) |> filter(n > 1) |> pull(Sample)
  if (length(dup) > 0) wca_abort("group_recurrence(): samples in more than one group: {str_c(dup, collapse = ', ')}")

  hits <- tbl |>
    select(Key = all_of(key_col), Sample, all_of(detail_cols)) |>
    inner_join(membership, by = "Sample")
  per_group <- hits |>
    distinct(Key, Group, Sample) |>
    group_by(Key, Group) |>
    summarize(n = n(), Samples = str_c(sort(Sample), collapse = ";"), .groups = "drop")
  wide <- purrr::reduce(names(groups), \(acc, g) {
    cols <- per_group |> filter(Group == g) |>
      transmute(Key, "n_{g}" := n, "PCT_{g}" := n / sizes[[g]], "Samples_{g}" := Samples)
    acc |> left_join(cols, by = "Key") |>
      mutate("n_{g}" := coalesce(.data[[str_c("n_", g)]], 0L),
        "PCT_{g}" := coalesce(.data[[str_c("PCT_", g)]], 0),
        "Samples_{g}" := coalesce(.data[[str_c("Samples_", g)]], ""))
  }, .init = hits |> distinct(Key, Sample) |> count(Key, name = "n"))
  if (length(detail_cols) > 0) {
    details <- hits |>
      group_by(Key) |>
      summarize(across(all_of(detail_cols), \(x) str_c(sort(unique(x)), collapse = ";")), .groups = "drop")
    wide <- left_join(wide, details, by = "Key")
  }
  out <- wide |>
    arrange(desc(n), across(all_of(str_c("n_", names(groups))), desc), Key) |>
    rename(!!key_col := Key)
  attr(out, "denominator") <- sizes
  out
}

#' Data dictionary rows for `group_recurrence()` columns
#'
#' @param groups The `groups` list given to `group_recurrence()`.
#' @return Tibble: Sheet, Column, Description, to bind onto
#'   `wca_dictionary()` before `write_report_xlsx()`.
group_recurrence_dictionary <- function(groups) {
  purrr::imap(groups, \(x, g) tribble(
    ~Sheet, ~Column, ~Description,
    "*", str_c("n_", g), str_glue("{g} samples with at least one event for the row; each sample counts once"),
    "*", str_c("PCT_", g), str_glue("n_{g} / number of {g} samples ({length(unique(x))})"),
    "*", str_c("Samples_", g), str_glue("{g} samples carrying the event")
  )) |> bind_rows() |> mutate(Description = as.character(Description))
}
