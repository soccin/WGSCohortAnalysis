# 20_genome.R
#
# Genome reference derived entirely from the bundled UCSC cytoBand table:
# chromosome order and lengths, centromeres (acen bands), arm boundaries
# and the acrocentric list. Every annotator takes `genome =` so a second
# build only needs an entry in `wca_genomes` and a cytoBand file.

wca_genomes <- list(
  hg19 = list(
    cytoband = "cytoBand_hg19.txt.gz",
    build_aliases = c("GRCh37", "hg19", "b37"),
    acrocentric = c("13", "14", "15", "21", "22")
  )
)

.wca_genome_cache <- new.env(parent = emptyenv())

#' Normalize chromosome names
#'
#' Strips a `chr` prefix and maps FACETS numeric codes: 23 to X, 24 to Y,
#' 25 and M to MT. Everything else passes through as character.
#'
#' @param x Character or numeric vector.
#' @return Character vector.
normalize_chrom <- function(x) {
  x <- str_remove(as.character(x), "^chr")
  case_when(
    x == "23" ~ "X",
    x == "24" ~ "Y",
    x == "25" ~ "MT",
    x == "M" ~ "MT",
    TRUE ~ x
  )
}

#' Chromosome levels in karyotype order
#'
#' @param genome Genome name.
#' @param include_mt Append MT.
#' @return Character vector `1..22, X, Y[, MT]`.
chrom_levels <- function(genome = "hg19", include_mt = TRUE) {
  lv <- genome_info(genome)$chrom_levels
  if (include_mt) c(lv, "MT") else lv
}

#' Chromosome factor in karyotype order
#'
#' @param x Chromosome names (normalized or not).
#' @param genome Genome name.
#' @return Factor with the genome's chromosome levels.
chrom_factor <- function(x, genome = "hg19") {
  factor(normalize_chrom(x), levels = chrom_levels(genome))
}

#' Autosome names
#'
#' @param genome Genome name.
#' @return Character vector `"1".."22"`.
autosomes <- function(genome = "hg19") {
  genome_info(genome)$autosomes
}

#' Keep autosomal rows only
#'
#' The one place the sex-chromosome rule lives. Every CNV function in the
#' toolkit calls this by default because FACETS does not model the tumor's
#' sex: a male X (one copy) is reported as a whole-chromosome loss and a
#' female X against a mis-set sex as a gain, which floods the recurrence
#' tables with events that are not somatic. Y is never emitted by FACETS
#' but is dropped here as well. Chromosome names are normalized first, so
#' the FACETS code 23 is caught.
#'
#' @param df Tibble with a chromosome column.
#' @param chrom Name of that column.
#' @param genome Genome name.
#' @return `df` restricted to `autosomes(genome)`.
filter_autosomes <- function(df, chrom = "chrom", genome = "hg19") {
  require_cols(df, chrom, "filter_autosomes()")
  df |> filter(normalize_chrom(.data[[chrom]]) %in% autosomes(genome))
}

#' Genome reference tables
#'
#' Parsed once per session and cached. `cytobands` is BED-style (0-based,
#' half-open starts) with `Band` = chrom + band, e.g. `"8p23.1"`.
#'
#' @param genome Genome name; must be in `wca_genomes`.
#' @return List: name, build_aliases, chrom_levels, autosomes, chrom_len
#'   (named numeric), cytobands, centromeres, arms, acrocentric.
genome_info <- function(genome = "hg19") {
  if (!genome %in% names(wca_genomes)) {
    wca_abort("Unknown genome '{genome}'; known: {str_c(names(wca_genomes), collapse = ', ')}")
  }
  if (!is.null(.wca_genome_cache[[genome]])) return(.wca_genome_cache[[genome]])
  spec <- wca_genomes[[genome]]

  cyto <- read_tsv(wca_file("data", spec$cytoband),
    col_names = c("chrom", "bandStart", "bandEnd", "band", "stain"),
    col_types = "cddcc", progress = FALSE) |>
    mutate(chrom = normalize_chrom(chrom)) |>
    filter(!str_detect(chrom, "_")) |>
    mutate(Band = str_c(chrom, band), arm_side = str_sub(band, 1, 1))

  levels <- c(as.character(1:22), "X", "Y")
  levels <- levels[levels %in% cyto$chrom]
  cyto <- cyto |>
    mutate(chrom = factor(chrom, levels = levels)) |>
    arrange(chrom, bandStart) |>
    mutate(chrom = as.character(chrom))

  chrom_len <- cyto |>
    group_by(chrom) |>
    summarize(len = max(bandEnd), .groups = "drop")
  chrom_len <- set_names(chrom_len$len, chrom_len$chrom)[levels]

  centromeres <- cyto |>
    filter(stain == "acen") |>
    group_by(chrom) |>
    summarize(cen_start = min(bandStart), cen_end = max(bandEnd), .groups = "drop") |>
    arrange(factor(chrom, levels = levels))

  arms <- cyto |>
    group_by(chrom, side = arm_side) |>
    summarize(start = min(bandStart), end = max(bandEnd), .groups = "drop") |>
    mutate(arm = str_c(chrom, side), length = end - start) |>
    filter(!(chrom %in% spec$acrocentric & side == "p")) |>
    mutate(chrom = factor(chrom, levels = levels)) |>
    arrange(chrom, side) |>
    mutate(chrom = as.character(chrom)) |>
    select(arm, chrom, side, start, end, length)

  info <- list(
    name = genome,
    build_aliases = spec$build_aliases,
    chrom_levels = levels,
    autosomes = as.character(1:22),
    chrom_len = chrom_len,
    cytobands = select(cyto, chrom, bandStart, bandEnd, band, stain, Band),
    centromeres = centromeres,
    arms = arms,
    acrocentric = spec$acrocentric
  )
  assign(genome, info, envir = .wca_genome_cache)
  info
}

#' Chromosome arm table
#'
#' @param genome Genome name.
#' @return Tibble: arm, chrom, side, start, end, length. No p arm for the
#'   acrocentric chromosomes.
chrom_arms <- function(genome = "hg19") {
  genome_info(genome)$arms
}

#' Annotate positions with their cytoband
#'
#' Joins on chromosome and `pos` in `(bandStart, bandEnd]`.
#'
#' @param df Data frame.
#' @param chrom,pos Column names (strings) holding the chromosome and the
#'   position.
#' @param genome Genome name.
#' @param name Name of the added band column (e.g. `"Band"`, `"bandA"`).
#' @param keep_bounds Also add `<name>_start` and `<name>_end`.
#' @return `df` with the band column(s) added; positions outside every band
#'   get `NA`.
annotate_cytoband <- function(df, chrom = "chrom", pos = "gene_start", genome = "hg19",
                              name = "Band", keep_bounds = FALSE) {
  bands <- genome_info(genome)$cytobands |>
    select(.chrom = chrom, .bstart = bandStart, .bend = bandEnd, .Band = Band)
  out <- df |>
    mutate(.chrom = normalize_chrom(.data[[chrom]]), .pos = as.numeric(.data[[pos]])) |>
    left_join(bands, by = join_by(.chrom, between(.pos, .bstart, .bend, bounds = "(]"))) |>
    select(-.chrom, -.pos)
  out[[name]] <- out$.Band
  if (keep_bounds) {
    out[[str_c(name, "_start")]] <- out$.bstart
    out[[str_c(name, "_end")]] <- out$.bend
  }
  select(out, -.Band, -.bstart, -.bend)
}

#' Chromosome arm of positions
#'
#' @inheritParams annotate_cytoband
#' @param name Name of the added arm column.
#' @return `df` with an arm column such as `"17p"`.
annotate_arm <- function(df, chrom = "chrom", pos = "gene_start", genome = "hg19", name = "arm") {
  arms <- chrom_arms(genome) |> select(.chrom = chrom, .astart = start, .aend = end, .arm = arm)
  out <- df |>
    mutate(.chrom = normalize_chrom(.data[[chrom]]), .pos = as.numeric(.data[[pos]])) |>
    left_join(arms, by = join_by(.chrom, between(.pos, .astart, .aend, bounds = "(]"))) |>
    select(-.chrom, -.pos, -.astart, -.aend)
  out[[name]] <- out$.arm
  select(out, -.arm)
}
