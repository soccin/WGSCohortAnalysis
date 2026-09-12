# 08_read_mei.R
#
# Polymorphic mobile element insertions, read from a population SV catalogue
# such as the 1000 Genomes phase 3 integrated SV map.
#
# These matter to an SV report for one reason. The element is absent from the
# reference, so reads that run off the end of it have to be placed somewhere
# else, and an aligner puts them on whatever the reference offers that looks
# closest. For an SVA or a LINE1 that is a poly-A tract, of which the genome
# has thousands. The mate stays behind, the pair ends up megabases apart
# pointing the way a deletion requires, and every caller reports a
# rearrangement. Intersecting called breakends with a catalogue of these
# insertions is the cheapest way to find out how many calls are that.

#' Read mobile element insertions from a population SV VCF
#'
#' Site records only; the genotype columns are dropped, so a multi-gigabyte
#' catalogue with thousands of samples still reads in seconds.
#'
#' @param vcf Path to the VCF, optionally gzipped.
#' @param types Element classes to keep. The default is the three insertion
#'   classes; `NULL` keeps every record.
#' @return Tibble: chrom, pos, end, id, svtype, af, ac, an, svlen, tsd, alt.
#'   `end` comes from the INFO `END` field; insertion records have no span
#'   and leave it equal to `pos`, but deletion and duplication records in the
#'   same file do carry one, which is how a record's size is recovered when
#'   `SVLEN` is absent.
#'   `tsd` is the target site duplication when the catalogue reports one,
#'   which is the sequence a read crossing the insertion is clipped at.
#'
#'   Catalogues disagree about where the element class lives. The 1000
#'   Genomes map puts it in `SVTYPE` (`SVTYPE=SVA`); gnomAD puts every
#'   insertion under `SVTYPE=INS` and names the class in the ALT allele
#'   (`<INS:ME:SVA>`). `svtype` is normalised from whichever carries it, so
#'   the two read the same way.
read_mei_catalogue <- function(vcf, types = c("ALU", "LINE1", "SVA")) {
  if (!fs::file_exists(vcf)) wca_abort("MEI catalogue not found: {vcf}")
  cmd <- str_glue("zcat -f {shQuote(vcf)} | grep -v '^#' | cut -f1,2,3,5,8")
  lines <- system(cmd, intern = TRUE)
  if (length(lines) == 0) wca_abort("No records read from {vcf}")
  info_field <- function(x, key) {
    str_extract(x, str_glue("(^|;){key}=([^;]*)"), group = 2)
  }
  # A few catalogue records carry a comma-separated value per alt allele.
  # Take the first rather than letting as.numeric() coerce to NA with a warning.
  info_num <- function(x, key) {
    as.numeric(str_extract(info_field(x, key), "^[0-9.eE+-]+"))
  }
  out <- read_tsv(I(lines), col_names = c("chrom", "pos", "id", "alt", "info"),
                  col_types = cols(.default = "c"), progress = FALSE) |>
    mutate(
      chrom = normalize_chrom(chrom),
      pos = as.integer(pos),
      # gnomAD names the element class in the ALT allele, 1000 Genomes in
      # SVTYPE. Prefer the ALT when it carries one.
      svtype = coalesce(str_extract(alt, "<INS:ME:([A-Za-z0-9]+)>", group = 1),
                        info_field(info, "SVTYPE")),
      af = info_num(info, "AF"),
      ac = as.integer(info_num(info, "AC")),
      an = as.integer(info_num(info, "AN")),
      svlen = as.integer(info_num(info, "SVLEN")),
      end = coalesce(as.integer(info_num(info, "END")), pos),
      tsd = na_if(info_field(info, "TSD"), "null")
    ) |>
    select(chrom, pos, end, id, svtype, af, ac, an, svlen, tsd, alt)
  if (!is.null(types)) out <- out |> filter(svtype %in% types)
  wca_msg("  MEI catalogue: {nrow(out)} records, types {str_c(sort(unique(out$svtype)), collapse = ', ')}")
  out
}

#' Called breakends that sit on a catalogued insertion
#'
#' Both ends of every call are tested. A call is reported once per breakend
#' that lands near an element, so a junction joining two of them appears
#' twice.
#'
#' @param sv Tibble from `sv_events()`: UUID, Sample, sv_class, CHROM_A,
#'   START_A, CHROM_B, START_B.
#' @param mei From `read_mei_catalogue()`.
#' @param window Base pairs either side of the catalogued position.
#' @return Tibble: UUID, Sample, sv_class, end (`A`/`B`), breakend_chrom,
#'   breakend_pos, mei_id, svtype, af, tsd, mei_pos, distance_bp.
breakends_near_mei <- function(sv, mei, window = 500L) {
  require_cols(sv, c("UUID", "Sample", "CHROM_A", "START_A", "CHROM_B", "START_B"), "sv")
  require_cols(mei, c("chrom", "pos", "id", "svtype"), "mei")
  sv <- sv |> mutate(sv_class = if ("sv_class" %in% names(sv)) sv_class else NA_character_)
  ends <- bind_rows(
    sv |> transmute(UUID, Sample, sv_class, end = "A", chrom = CHROM_A, bp = START_A),
    sv |> transmute(UUID, Sample, sv_class, end = "B", chrom = CHROM_B, bp = START_B)
  ) |>
    filter(!is.na(chrom), !is.na(bp))
  # `mei$end` would collide with the A/B label column below and silently
  # resolve to base::end(), so drop it here: a breakend match only needs pos.
  ends |>
    inner_join(mei |> select(-any_of("end")), by = "chrom",
               relationship = "many-to-many") |>
    filter(abs(bp - pos) <= window) |>
    transmute(UUID, Sample, sv_class, end, breakend_chrom = chrom, breakend_pos = bp,
              mei_id = id, svtype, af, tsd, mei_pos = pos,
              distance_bp = as.integer(abs(bp - pos))) |>
    arrange(UUID, end, distance_bp)
}

#' Which of a set of positions fall inside a BED or BEDPE filter
#'
#' Answers whether a filter a pipeline already ships would have removed a
#' locus. A BEDPE is read as two independent interval sets, because a pair
#' filter only fires when both ends match and that is rarely what is wanted
#' when the question is about one locus.
#'
#' @param positions Tibble with `chrom` and `pos`.
#' @param file BED or BEDPE path, optionally gzipped.
#' @param slop Widen every interval by this many base pairs.
#' @return Logical vector, one per row of `positions`.
positions_in_filter <- function(positions, file, slop = 0L) {
  require_cols(positions, c("chrom", "pos"), "positions")
  if (!fs::file_exists(file)) wca_abort("Filter file not found: {file}")
  raw <- read_tsv(file, col_names = FALSE, col_types = cols(.default = "c"),
                  comment = "#", progress = FALSE)
  take <- function(cols) {
    raw |>
      select(all_of(cols)) |>
      set_names(c("chrom", "start", "end")) |>
      # A track line or a "." placeholder in a coordinate column is dropped
      # rather than coerced, so the filter never silently widens.
      mutate(chrom = normalize_chrom(chrom),
             start = suppressWarnings(as.numeric(start)) - slop,
             end = suppressWarnings(as.numeric(end)) + slop) |>
      filter(!is.na(start), !is.na(end))
  }
  iv <- if (ncol(raw) >= 6) bind_rows(take(1:3), take(4:6)) else take(1:3)
  by_chrom <- split(iv, iv$chrom)
  map2_lgl(positions$chrom, positions$pos, \(cc, pp) {
    h <- by_chrom[[cc]]
    !is.null(h) && any(h$start <= pp & h$end >= pp)
  })
}
