# 27_read_evidence.R
#
# What the alignments say, for the cases where a caller's summary is not
# enough to decide what happened. An SV caller reports a junction; the
# reads that produced it report where they actually sit, how far their
# mates went, how much of each read the aligner had to clip away, and what
# that clipped sequence is made of.
#
# Three questions this module answers directly:
#
#   Is a called deletion a real loss?  `bam_binned_counts()` over the
#   called interval and a flank either side. A heterozygous loss halves
#   the count; an artefact does not move it.
#
#   Where do the clips fall, and what are they?  `bam_soft_clips()` puts
#   one row per soft clip at its reference position with the clipped
#   sequence, so a cluster of clips at one base is visible as a cluster.
#
#   Is the locus a genotype?  `clip_genotype()` counts, per sample, reads
#   that cross the locus intact against reads clipped at it whose clipped
#   sequence matches a named class. Zero intact reads and a pile of clips
#   is a homozygous insertion; both present is a heterozygote.
#
# Everything here needs Rsamtools, which should not be attached with
# `library()`: it and its dependencies mask `count()` and friends, and the
# toolkit's dplyr pipelines break silently when they do. The functions
# below use `::`. CIGAR arithmetic is done here rather than through
# GenomicAlignments so the module does not pull in the Bioconductor
# sequence-info classes at all.

#' Are the alignment packages available
#'
#' @param what Name used in the message when they are not.
#' @return TRUE, or FALSE with a message.
have_bam_support <- function(what = "alignment evidence") {
  ok <- requireNamespace("Rsamtools", quietly = TRUE) &&
    requireNamespace("GenomicRanges", quietly = TRUE)
  if (!ok) wca_msg("  Rsamtools not installed; skipping {what}")
  ok
}

#' Reference bases a CIGAR consumes
#'
#' M, =, X, D and N advance along the reference; I, S, H and P do not.
#'
#' @param cigar Character vector of CIGAR strings.
#' @return Integer vector, one width per input.
cigar_ref_width <- function(cigar) {
  ops <- str_match_all(cigar, "([0-9]+)([MIDNSHP=X])")
  vapply(ops, function(m) {
    if (nrow(m) == 0) return(NA_integer_)
    sum(as.integer(m[, 2])[m[, 3] %in% c("M", "=", "X", "D", "N")])
  }, integer(1))
}

#' Default flag filter: primary, mapped, non-duplicate alignments
#'
#' @return A `scanBamFlag` object.
primary_reads_flag <- function() {
  Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                         isSecondaryAlignment = FALSE,
                         isSupplementaryAlignment = FALSE,
                         isDuplicate = FALSE)
}

#' Reads overlapping a window, as a tibble
#'
#' One row per primary alignment. `end` is the last reference base the
#' alignment covers, so `pos` and `end` bracket the aligned part only:
#' soft-clipped bases fall outside them and are reported by
#' `bam_soft_clips()`.
#'
#' @param bam Path to an indexed BAM.
#' @param chrom,start,end Window, one based and inclusive.
#' @param min_mapq Drop alignments below this mapping quality.
#' @return Tibble: read_id, qname, flag, chrom, pos, ref_end, strand, mapq,
#'   cigar, mate_chrom, mate_pos, isize, seq. Empty tibble when the window
#'   holds nothing.
bam_reads <- function(bam, chrom, start, end, min_mapq = 0L) {
  if (!have_bam_support("bam_reads()")) return(empty_bam_reads())
  if (!file_exists(bam)) wca_abort("BAM not found: {bam}")
  rng <- GenomicRanges::GRanges(as.character(chrom),
                                IRanges::IRanges(as.integer(start), as.integer(end)))
  param <- Rsamtools::ScanBamParam(
    which = rng,
    what = c("qname", "flag", "rname", "pos", "strand", "mapq", "cigar",
             "mrnm", "mpos", "isize", "seq"),
    flag = primary_reads_flag())
  b <- Rsamtools::scanBam(bam, param = param)[[1]]
  if (length(b$pos) == 0) return(empty_bam_reads())

  ref_width <- cigar_ref_width(b$cigar)
  out <- tibble(
    qname = b$qname,
    flag = as.integer(b$flag),
    chrom = as.character(b$rname),
    pos = as.integer(b$pos),
    ref_end = as.integer(b$pos) + ref_width - 1L,
    strand = as.character(b$strand),
    mapq = as.integer(b$mapq),
    cigar = b$cigar,
    mate_chrom = as.character(b$mrnm),
    mate_pos = as.integer(b$mpos),
    isize = as.integer(b$isize),
    seq = as.character(b$seq))
  out <- out |> filter(mapq >= min_mapq)
  out |> mutate(read_id = row_number(), .before = 1)
}

empty_bam_reads <- function() {
  tibble(read_id = integer(), qname = character(), flag = integer(),
         chrom = character(), pos = integer(), ref_end = integer(),
         strand = character(), mapq = integer(), cigar = character(),
         mate_chrom = character(), mate_pos = integer(), isize = integer(),
         seq = character())
}

#' Soft clips of a read table, one row per clip
#'
#' A read can be clipped at both ends, so this is long rather than wide.
#' `clip_pos` is the reference base the clip abuts: the first aligned base
#' for a leading clip, the last aligned base for a trailing one. Clipped
#' bases themselves have no reference position, which is the point.
#'
#' @param reads Output of `bam_reads()`.
#' @param min_len Ignore clips shorter than this.
#' @return Tibble: read_id, side (`left`/`right`), clip_pos, clip_len,
#'   clip_seq, gc_frac, at_frac, plus pos, ref_end, mapq carried through.
bam_soft_clips <- function(reads, min_len = 10L) {
  require_cols(reads, c("read_id", "cigar", "seq", "pos", "ref_end"), "reads")
  lead <- as.integer(str_match(reads$cigar, "^([0-9]+)S")[, 2])
  trail <- as.integer(str_match(reads$cigar, "([0-9]+)S$")[, 2])

  left <- reads |>
    mutate(side = "left", clip_len = lead, clip_pos = pos,
           clip_seq = str_sub(seq, 1L, coalesce(lead, 0L)))
  right <- reads |>
    mutate(side = "right", clip_len = trail, clip_pos = ref_end,
           clip_seq = str_sub(seq, nchar(seq) - coalesce(trail, 0L) + 1L, nchar(seq)))

  bind_rows(left, right) |>
    filter(!is.na(clip_len), clip_len >= min_len) |>
    mutate(gc_frac = str_count(clip_seq, "[GCgc]") / clip_len,
           at_frac = str_count(clip_seq, "[ATat]") / clip_len) |>
    select(read_id, side, clip_pos, clip_len, clip_seq, gc_frac, at_frac,
           pos, ref_end, mapq, strand, mate_chrom, mate_pos) |>
    arrange(clip_pos, side)
}

#' Read counts in bins across a region, for a copy-number sanity check
#'
#' Counts rather than depth: each bin is the same width, so the counts are
#' comparable without normalising for read length. Bins are sampled evenly
#' across the region rather than tiled, so a 50 Mb interval costs a few
#' seconds instead of an hour.
#'
#' @param bam Path to an indexed BAM.
#' @param chrom,start,end Region.
#' @param n_bins How many sampling windows.
#' @param bin_width Width of each window in bases.
#' @return Tibble: chrom, bin_start, bin_end, n, rel (n over the median n).
bam_binned_counts <- function(bam, chrom, start, end, n_bins = 120L, bin_width = 2000L) {
  if (!have_bam_support("bam_binned_counts()")) {
    return(tibble(chrom = character(), bin_start = integer(), bin_end = integer(),
                  n = integer(), rel = double()))
  }
  if (!file_exists(bam)) wca_abort("BAM not found: {bam}")
  starts <- as.integer(round(seq(start, end - bin_width, length.out = n_bins)))
  rng <- GenomicRanges::GRanges(as.character(chrom),
                                IRanges::IRanges(starts, width = bin_width))
  cb <- Rsamtools::countBam(bam, param = Rsamtools::ScanBamParam(
    which = rng, flag = primary_reads_flag()))
  n <- as.integer(cb$records)
  med <- stats::median(n)
  tibble(chrom = as.character(chrom), bin_start = starts,
         bin_end = starts + bin_width - 1L, n = n,
         rel = if (med > 0) n / med else NA_real_)
}

#' Per base depth across a small window
#'
#' Depth is walked off the CIGAR rather than taken from a GAlignments
#' object, so deletions and skips are not counted as covered and the
#' function does not depend on the Bioconductor sequence-info classes,
#' whose S4 dispatch is fragile across installs.
#'
#' @param bam Path to an indexed BAM.
#' @param chrom,start,end Window; keep it under a few tens of kilobases.
#' @param min_mapq Drop alignments below this mapping quality.
#' @return Tibble: chrom, pos, depth.
bam_base_depth <- function(bam, chrom, start, end, min_mapq = 0L) {
  start <- as.integer(start)
  end <- as.integer(end)
  reads <- bam_reads(bam, chrom, start, end, min_mapq = min_mapq)
  width <- end - start + 1L
  delta <- integer(width + 2L)
  if (nrow(reads) > 0) {
    ops <- str_match_all(reads$cigar, "([0-9]+)([MIDNSHP=X])")
    for (i in seq_len(nrow(reads))) {
      ref <- reads$pos[i]
      m <- ops[[i]]
      for (j in seq_len(nrow(m))) {
        n <- as.integer(m[j, 2])
        op <- m[j, 3]
        if (op %in% c("M", "=", "X")) {
          a <- max(ref, start)
          b <- min(ref + n - 1L, end)
          if (b >= a) {
            delta[a - start + 1L] <- delta[a - start + 1L] + 1L
            delta[b - start + 2L] <- delta[b - start + 2L] - 1L
          }
          ref <- ref + n
        } else if (op %in% c("D", "N")) {
          ref <- ref + n
        }
      }
    }
  }
  tibble(chrom = as.character(chrom),
         pos = seq.int(start, end),
         depth = cumsum(delta)[seq_len(width)])
}

#' Genotype a locus from intact reads against clipped reads
#'
#' At an insertion site, a read either crosses the site in the reference
#' (the allele without the insert) or runs into sequence that is not in
#' the reference and gets clipped (the allele with it). Counting the two
#' is a genotype: no intact reads means homozygous, both means
#' heterozygous, no clips means the reference allele.
#'
#' @param bam Path to an indexed BAM.
#' @param chrom Chromosome.
#' @param start,end The locus, usually the target site duplication. An
#'   intact read must cover both ends plus `flank`.
#' @param classes Named list of predicates on the clipped sequence, each
#'   taking a character vector and returning logical. One count column per
#'   name.
#' @param window Bases either side of the locus to fetch.
#' @param min_clip Ignore clips shorter than this.
#' @param flank Bases beyond the locus an intact read must also cover.
#' @return One row tibble: n_reads, n_intact, n_clipped, one column per
#'   class, and vaf_clipped = n_clipped / (n_clipped + n_intact).
clip_genotype <- function(bam, chrom, start, end, classes = list(),
                          window = 300L, min_clip = 15L, flank = 10L) {
  reads <- bam_reads(bam, chrom, start - window, end + window)
  base <- tibble(n_reads = nrow(reads), n_intact = 0L, n_clipped = 0L)
  if (nrow(reads) == 0) {
    for (nm in names(classes)) base[[nm]] <- 0L
    return(base |> mutate(vaf_clipped = NA_real_))
  }
  intact <- reads |>
    filter(!str_detect(cigar, "S"),
           pos <= start - flank, ref_end >= end + flank)
  clips <- bam_soft_clips(reads, min_len = min_clip) |>
    filter(clip_pos >= start - flank, clip_pos <= end + flank)

  out <- tibble(n_reads = nrow(reads), n_intact = nrow(intact),
                n_clipped = nrow(clips))
  for (nm in names(classes)) {
    hit <- if (nrow(clips) == 0) logical(0) else classes[[nm]](clips$clip_seq)
    out[[nm]] <- sum(hit, na.rm = TRUE)
  }
  out |> mutate(vaf_clipped = if (n_clipped + n_intact > 0)
    n_clipped / (n_clipped + n_intact) else NA_real_)
}

#' Fraction of A and T in sliding reference windows
#'
#' A breakend that sits in a trough of G+C is a breakend the aligner
#' placed on composition rather than on sequence identity.
#'
#' @param fasta Indexed reference fasta.
#' @param chrom,start,end Region.
#' @param window Window width in bases.
#' @param step Step between windows.
#' @return Tibble: chrom, pos (window midpoint), at_frac.
ref_at_profile <- function(fasta, chrom, start, end, window = 50L, step = 10L) {
  seqx <- read_ref_seq(fasta, chrom, start, end)
  if (is.na(seqx)) return(tibble(chrom = character(), pos = integer(), at_frac = double()))
  offs <- seq.int(1L, nchar(seqx) - window + 1L, by = step)
  win <- str_sub(seqx, offs, offs + window - 1L)
  tibble(chrom = as.character(chrom),
         pos = as.integer(start) + offs - 1L + window %/% 2L,
         at_frac = str_count(win, "[ATat]") / window)
}

#' Read pair orientation, from the SAM flag
#'
#' The orientation of a discordant pair is what separates the classes of
#' rearrangement, so it is worth reading off the flag rather than trusting
#' a caller's label. Orientation is expressed from the leftmost mate:
#'
#'   FR  the normal orientation. A deletion produces it, just with a
#'       larger insert size than the library.
#'   RF  the orientation a tandem duplication produces.
#'   FF  one of the two junctions of an inversion (head to head).
#'   RR  the other junction of an inversion (tail to tail).
#'
#' A balanced inversion produces FF and RR pairs and no FR pairs across
#' its span. Finding only FR pairs between two loci rules an inversion of
#' the segment between them out, whatever else may be going on.
#'
#' @param reads Output of `bam_reads()`; needs flag, pos and mate_pos.
#' @return `reads` with self_reverse, mate_reverse, self_is_left,
#'   orientation (`FR`/`RF`/`FF`/`RR`) and sv_signature, the class of
#'   rearrangement that orientation is consistent with.
read_pair_orientation <- function(reads) {
  require_cols(reads, c("flag", "pos", "mate_pos"), "reads")
  self_rev <- bitwAnd(reads$flag, 16L) > 0
  mate_rev <- bitwAnd(reads$flag, 32L) > 0
  self_left <- !is.na(reads$mate_pos) & reads$pos <= reads$mate_pos
  left_rev <- if_else(self_left, self_rev, mate_rev)
  right_rev <- if_else(self_left, mate_rev, self_rev)
  orient <- str_c(if_else(left_rev, "R", "F"), if_else(right_rev, "R", "F"))
  orient[is.na(reads$mate_pos)] <- NA_character_
  reads |>
    mutate(self_reverse = self_rev, mate_reverse = mate_rev,
           self_is_left = self_left,
           orientation = orient,
           sv_signature = recode(orient,
             FR = "deletion or normal", RF = "tandem duplication",
             FF = "inversion, head to head", RR = "inversion, tail to tail",
             .default = NA_character_))
}

#' How much of a clipped sequence actually comes from a candidate locus
#'
#' A junction between two loci requires that the sequence clipped off at
#' one of them is the sequence found at the other. Which orientation it
#' should be in is the thing under test: a deletion joins the partner as
#' the reference reads it, an inversion joins its reverse complement. Both
#' are checked, so a clipped sequence that matches neither belongs to
#' neither model.
#'
#' @param clips Output of `bam_soft_clips()`.
#' @param partner_seq Reference sequence at the candidate partner locus.
#' @return `clips` with match_forward, match_revcomp (longest exact match
#'   in bases) and match_frac, the better of the two over the clip length.
clip_partner_match <- function(clips, partner_seq) {
  require_cols(clips, c("clip_seq", "clip_len"), "clips")
  if (is.na(partner_seq) || nchar(partner_seq) == 0) {
    wca_abort("clip_partner_match(): partner_seq is empty")
  }
  rc <- revcomp(partner_seq)
  clips |>
    mutate(
      match_forward = purrr::map_int(clip_seq, \(s) longest_common_substring(s, partner_seq)$len),
      match_revcomp = purrr::map_int(clip_seq, \(s) longest_common_substring(s, rc)$len),
      match_frac = pmax(match_forward, match_revcomp) / clip_len)
}
