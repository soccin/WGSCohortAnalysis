# 03_read_maf.R
#
# Tempo somatic MAF reader. Columns are selected by name so the 259 and 270
# column schemas (the newer one adds neo_* columns) read identically. Rows
# are filtered per file before binding so a 200-pair cohort fits in memory.

#' Standard MAF columns kept by `read_tempo_maf()`
#'
#' @return Character vector read from `data/maf_standard_cols.txt`.
maf_standard_cols <- function() {
  read_list_file(wca_file("data", "maf_standard_cols.txt"))
}

#' Variant classes that are neither coding nor splice
#'
#' @return Character vector.
maf_noncoding_classes <- function() {
  c("Intron", "IGR", "3'Flank", "3'UTR", "5'Flank", "5'UTR", "Silent", "RNA")
}

#' Loss-of-function variant classes
#'
#' @return Character vector.
maf_lof_classes <- function() {
  c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site",
    "Translation_Start_Site", "Nonstop_Mutation")
}

#' The canonical non-silent row filter
#'
#' Drops silent and non-coding classes but keeps 5'Flank rows for the genes
#' in `keep_5p_flank` (the TERT promoter hotspot has no HGVSp annotation and
#' would otherwise be lost). Returns a function so it can be passed as
#' `row_filter` and cached by content.
#'
#' @param keep_5p_flank Genes whose 5'Flank rows are retained.
#' @return A function `maf -> maf`.
non_silent <- function(keep_5p_flank = "TERT") {
  bad <- maf_noncoding_classes()
  keep <- keep_5p_flank
  structure(
    function(maf) {
      filter(maf,
        (Hugo_Symbol %in% keep & Variant_Classification == "5'Flank") |
          !Variant_Classification %in% bad)
    },
    wca_filter = list(name = "non_silent", keep_5p_flank = keep)
  )
}

#' Row filter that keeps every row
#'
#' @return A function `maf -> maf`.
keep_all_rows <- function() {
  structure(function(maf) maf, wca_filter = list(name = "keep_all_rows"))
}

#' Read one Tempo somatic MAF
#'
#' @param file Path to `*.somatic.final.maf` (plain or gzipped).
#' @param cols Columns to keep, selected by name; missing ones are warned
#'   about, not fatal. `NULL` keeps every column.
#' @param row_filter Function applied to the tibble after reading; see
#'   `non_silent()`. `NULL` keeps all rows.
#' @param genome Expected build; `NCBI_Build` must match one of its aliases.
#' @param add_uuid Add `UUID` = Sample:Chrom:Start:Ref:Alt.
#' @return Typed tibble.
read_tempo_maf <- function(file, cols = maf_standard_cols(), row_filter = non_silent(),
                           genome = "hg19", add_uuid = TRUE) {
  n_skip <- sum(str_starts(readr::read_lines(file, n_max = 5), "#"))
  header <- vroom::vroom(file, delim = "\t", n_max = 0, skip = n_skip,
    col_types = cols(.default = "c"), show_col_types = FALSE, progress = FALSE)
  have <- names(header)
  if (!is.null(cols)) {
    missing <- setdiff(cols, have)
    if (length(missing) > 0) {
      wca_warn("{basename(file)}: {length(missing)} requested MAF columns absent: {str_c(missing, collapse = ', ')}")
    }
    keep <- intersect(cols, have)
  } else {
    keep <- have
  }
  maf <- vroom::vroom(file, delim = "\t", col_select = all_of(keep), skip = n_skip,
    col_types = cols(.default = "c"), na = c("", "NA", "."), quote = "",
    show_col_types = FALSE, progress = FALSE)
  if (!is.null(row_filter)) maf <- row_filter(maf)
  if (nrow(maf) == 0) return(type_convert_silent(maf))
  if ("NCBI_Build" %in% names(maf)) {
    builds <- unique(maf$NCBI_Build)
    ok <- genome_info(genome)$build_aliases
    if (!all(builds %in% ok)) {
      wca_abort("{basename(file)}: NCBI_Build {str_c(builds, collapse = ',')} is not {genome} ({str_c(ok, collapse = '/')})")
    }
  }
  maf <- type_convert_silent(maf)
  if ("Chromosome" %in% names(maf)) maf <- mutate(maf, Chromosome = normalize_chrom(Chromosome))
  if (add_uuid) {
    maf <- mutate(maf, UUID = str_c(Tumor_Sample_Barcode, Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2, sep = ":"))
  }
  maf
}

#' Read the MAF of every cohort pair
#'
#' @param cohort Tibble from `build_cohort()`.
#' @param cache_dir Directory for `cached_read()`; `NULL` disables caching.
#' @param ... Passed to `read_tempo_maf()`.
#' @return Bound tibble with `Sample` and `NID` from the cohort first.
read_cohort_maf <- function(cohort, cache_dir = NULL, ...) {
  read_cohort_files(cohort, file_col = "snv_file", reader = read_tempo_maf,
    cache_dir = cache_dir, what = "MAF", ...)
}

#' Map a reader over cohort files and bind the results
#'
#' @param cohort Cohort tibble.
#' @param file_col Column of `cohort` holding the file paths, or a named
#'   character vector of paths keyed by `Sample`.
#' @param reader Reader function taking the path first.
#' @param cache_dir Passed to `cached_read()`; `NULL` reads directly.
#' @param what Label for progress messages.
#' @param ... Passed to the reader.
#' @return Bound tibble with `Sample`, `NID` first; empty tibble when no
#'   file yields rows.
read_cohort_files <- function(cohort, file_col, reader, cache_dir = NULL, what = "file", ...) {
  reader_name <- deparse1(substitute(reader))
  if (is.character(file_col) && length(file_col) == 1 && file_col %in% names(cohort)) {
    files <- cohort |> select(Sample, NID, file = all_of(file_col))
  } else {
    files <- tibble(Sample = names(file_col), file = unname(file_col)) |>
      left_join(select(cohort, Sample, NID), by = "Sample")
  }
  files <- filter(files, !is.na(file))
  exists <- fs::file_exists(files$file)
  wca_msg("  {what}: {sum(exists)} files found, {sum(!exists)} missing, {nrow(cohort) - nrow(files)} pairs without a path")
  files <- files[exists, ]
  if (nrow(files) == 0) return(tibble(Sample = character(), NID = character()))
  read1 <- function(f) {
    if (is.null(cache_dir)) reader(f, ...) else
      cached_read(f, reader, ..., cache_dir = cache_dir, reader_name = reader_name)
  }
  tabs <- purrr::map(files$file, read1, .progress = !isTRUE(getOption("wca.quiet", FALSE)))
  names(tabs) <- files$Sample
  tabs <- purrr::compact(tabs)
  tabs <- tabs[purrr::map_int(tabs, nrow) > 0]
  if (length(tabs) == 0) return(tibble(Sample = character(), NID = character()))
  # ids parsed from the file itself are kept as *_file; the cohort ids win
  tabs <- purrr::map(tabs, \(x) rename(x, any_of(c(Sample_file = "Sample", NID_file = "NID"))))
  out <- data.table::rbindlist(tabs, fill = TRUE, idcol = "Sample") |>
    as_tibble() |>
    left_join(select(files, Sample, NID), by = "Sample") |>
    select(Sample, NID, everything())
  n_rows <- nrow(out)
  wca_msg("  {what}: {n_rows} rows across {n_distinct(out$Sample)} samples")
  out
}
