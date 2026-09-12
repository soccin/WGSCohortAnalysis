# 26_junction_seq.R
#
# What a called SV junction is actually made of. An SV caller reports two
# breakends; the assembled consensus that spans them is the evidence. If
# each half of that consensus matches unique sequence at the locus it was
# assigned to, the junction is real. If one half is a homopolymer, that
# breakend carries no information and the caller has parked it wherever
# the aligner scored best.
#
# `junction_anchor_support()` is the entry point. Reference sequence needs
# Rsamtools and Biostrings; everything else is plain string work, so the
# triage runs without Bioconductor if the sequence is supplied directly.

#' Longest exact common substring of two sequences
#'
#' Binary search on the length, so a 200 bp consensus against a 1 kb
#' window costs a handful of passes rather than the naive product.
#'
#' @param a,b Character scalars.
#' @return List: len, a_start, b_start, seq. `len` is 0 when nothing matches.
longest_common_substring <- function(a, b) {
  none <- list(len = 0L, a_start = NA_integer_, b_start = NA_integer_, seq = "")
  if (is.na(a) || is.na(b) || nchar(a) == 0 || nchar(b) == 0) return(none)
  find_of_length <- function(k) {
    for (i in seq_len(nchar(a) - k + 1)) {
      s <- str_sub(a, i, i + k - 1)
      hit <- str_locate(b, fixed(s))
      if (!is.na(hit[1, 1])) {
        return(list(len = k, a_start = i, b_start = unname(hit[1, 1]), seq = s))
      }
    }
    NULL
  }
  best <- none
  lo <- 1L
  hi <- nchar(a)
  while (lo <= hi) {
    mid <- (lo + hi) %/% 2L
    got <- find_of_length(mid)
    if (is.null(got)) hi <- mid - 1L else { best <- got; lo <- mid + 1L }
  }
  best
}

#' Reverse complement of a DNA string
#'
#' Needed whenever a junction has to be tested in both orientations: a
#' deletion joins two loci as the reference reads them, an inversion joins
#' one of them reverse complemented, and only the sequence can tell the
#' two apart.
#'
#' @param x Character vector of sequences.
#' @return Character vector, same length. Non-ACGT characters are left as
#'   they are apart from case-matched complementing.
revcomp <- function(x) {
  stringi::stri_reverse(chartr("ACGTacgtNn", "TGCAtgcaNn", x))
}

#' Reference sequence for a window, from an indexed fasta
#'
#' @param fasta Path to a fasta with a `.fai` index.
#' @param chrom,start,end Window, one based and inclusive.
#' @return Character scalar, or NA if Rsamtools is not installed.
read_ref_seq <- function(fasta, chrom, start, end) {
  if (!requireNamespace("Rsamtools", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE)) {
    wca_msg("  Rsamtools not installed; skipping reference sequence")
    return(NA_character_)
  }
  if (!file_exists(fasta)) wca_abort("Reference fasta not found: {fasta}")
  rng <- GenomicRanges::GRanges(as.character(chrom),
                                IRanges::IRanges(as.integer(start), as.integer(end)))
  as.character(Rsamtools::scanFa(Rsamtools::FaFile(fasta), rng))
}

#' How informative each end of a junction consensus is
#'
#' Matches the consensus against a reference window for each breakend and
#' reports, per breakend, the longest exact match and how much of it is
#' anything other than A and T. A poly-A or (TAAAA)n tract has a GC count
#' of zero and occurs in thousands of places, so a breakend anchored only
#' on one is not located by the data.
#'
#' @param consensus Character scalar, the junction-spanning contig.
#' @param loci Tibble with breakend, chrom, pos, and either `seq` or
#'   enough for `read_ref_seq()` given `fasta`.
#' @param fasta Optional indexed fasta used when `loci$seq` is absent.
#' @param window Bases either side of `pos` to pull.
#' @param min_gc Informative anchors need at least this many G or C bases.
#' @return Tibble: breakend, chrom, pos, match_len, match_start_in_consensus,
#'   match_pos, gc_in_match, at_fraction, informative.
junction_anchor_support <- function(consensus, loci, fasta = NULL, window = 300L,
                                    min_gc = 5L) {
  require_cols(loci, c("breakend", "chrom", "pos"), "loci")
  if (!"seq" %in% names(loci)) {
    if (is.null(fasta)) wca_abort("junction_anchor_support() needs loci$seq or fasta =")
    loci <- loci |>
      mutate(seq = pmap_chr(list(chrom, pos), \(cc, pp) {
        read_ref_seq(fasta, cc, pp - window, pp + window)
      }))
  }
  loci |>
    mutate(hit = map(seq, \(s) longest_common_substring(consensus, s))) |>
    mutate(
      match_len = map_int(hit, \(h) as.integer(h$len)),
      match_start_in_consensus = map_int(hit, \(h) as.integer(h$a_start)),
      match_pos = as.integer(pos - window + map_int(hit, \(h) as.integer(h$b_start)) - 1L),
      matched_seq = map_chr(hit, \(h) h$seq),
      gc_in_match = str_count(matched_seq, "[GCgc]"),
      at_fraction = if_else(match_len > 0, 1 - gc_in_match / pmax(match_len, 1), NA_real_),
      informative = match_len > 0 & gc_in_match >= min_gc
    ) |>
    select(breakend, chrom, pos, match_len, match_start_in_consensus, match_pos,
           gc_in_match, at_fraction, informative, matched_seq)
}

#' Draw where each breakend's evidence sits along the consensus
#'
#' @param anchors Result of `junction_anchor_support()`, one set per sample
#'   if a `Sample` column is present.
#' @param consensus_len Named integer vector of consensus lengths, or a scalar.
#' @param title Plot title.
#' @return A ggplot object.
plot_junction_anatomy <- function(anchors, consensus_len, title = NULL) {
  # `anchors` usually carries its own consensus_len column, which would
  # mask the argument inside mutate(), so resolve the lengths first.
  lens <- consensus_len
  dat <- anchors |>
    mutate(Sample = if ("Sample" %in% names(anchors)) Sample else "consensus") |>
    mutate(len = if (length(lens) == 1) unname(lens) else unname(lens[Sample]),
           xmin = match_start_in_consensus,
           xmax = match_start_in_consensus + match_len - 1L,
           label = str_glue("{chrom}:{pos} ({gc_in_match} GC)"))
  if (anyNA(dat$len)) {
    wca_abort("plot_junction_anatomy(): no consensus length for {str_c(unique(dat$Sample[is.na(dat$len)]), collapse = ', ')}")
  }
  ggplot(dat) +
    geom_segment(aes(x = 1, xend = len, y = Sample, yend = Sample),
                 colour = "grey80", linewidth = 3) +
    geom_segment(aes(x = xmin, xend = xmax, y = Sample, yend = Sample,
                     colour = informative), linewidth = 5) +
    geom_text(aes(x = (xmin + xmax) / 2, y = Sample, label = gc_in_match),
              size = 2.4, colour = "white") +
    scale_colour_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d95f02"),
                        labels = c(`TRUE` = "unique sequence",
                                   `FALSE` = "A/T tract, not located"),
                        name = NULL) +
    labs(title = title, subtitle = "numbers are G or C bases in the matched segment",
         x = "position along the assembled junction consensus", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())
}
