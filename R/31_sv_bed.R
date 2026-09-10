# 31_sv_bed.R
#
# IGV breakpoint BED from an SV table (from dumpSVBED.R). Each SV yields two
# single-base features, one per breakend, named Sample_class_GenePair.

#' Breakpoint BED table from an SV table
#'
#' @param sv Tibble from `sv_events()` or `read_tempo_sv()` with
#'   `Sample` (or `TUMOR_ID`), `GenePair`, `CHROM_A/B`, `START_A/B`,
#'   `END_A/B`, `STRAND_A/B` (or `STRANDS`) and `sv_class`.
#' @param score BED score column value.
#' @param color BED itemRgb value.
#' @return Tibble in BED9 order: chrom, start, end, name, score, strand,
#'   thickStart, thickEnd, itemRgb.
sv_to_bed <- function(sv, score = 1000, color = "0,0,0") {
  if (!"Sample" %in% names(sv)) sv <- mutate(sv, Sample = TUMOR_ID)
  if (!"STRAND_A" %in% names(sv)) {
    sv <- mutate(sv, STRAND_A = str_sub(STRANDS, 1, 1), STRAND_B = str_sub(STRANDS, 2, 2))
  }
  if (!"sv_class" %in% names(sv)) sv <- mutate(sv, sv_class = sv_class(ID))
  base <- sv |>
    transmute(name = str_c(Sample, sv_class, GenePair, sep = "_"),
      CHROM_A, START_A, END_A, STRAND_A, CHROM_B, START_B, END_B, STRAND_B)
  a <- base |> transmute(chrom = CHROM_A, start = START_A, end = END_A, name, strand = STRAND_A)
  b <- base |> transmute(chrom = CHROM_B, start = START_B, end = END_B, name, strand = STRAND_B)
  bind_rows(a, b) |>
    mutate(score = score, thickStart = start, thickEnd = start, itemRgb = color) |>
    select(chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb) |>
    arrange(name, chrom, start)
}

#' Write a breakpoint BED file
#'
#' @param sv SV table (see `sv_to_bed()`).
#' @param file Output path.
#' @param track_name IGV track name written in the header line.
#' @return `file`, invisibly.
write_sv_bed <- function(sv, file, track_name = "SV breakpoints") {
  bed <- sv_to_bed(sv)
  readr::write_lines(str_glue('track name="{track_name}" itemRgb="On"'), file)
  readr::write_tsv(bed, file, col_names = FALSE, append = TRUE)
  invisible(file)
}
