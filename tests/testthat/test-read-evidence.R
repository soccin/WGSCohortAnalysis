# A synthetic BAM, built from a SAM written here, so the tests own their
# input and no sequencing data of any kind lives in the repository.

skip_without_bam <- function() {
  testthat::skip_if_not_installed("Rsamtools")
}

test_that("cigar_ref_width counts only the operations that consume reference", {
  expect_equal(cigar_ref_width("50M"), 50L)
  expect_equal(cigar_ref_width("10S40M"), 40L)
  expect_equal(cigar_ref_width("40M10S"), 40L)
  expect_equal(cigar_ref_width("20M5D20M"), 45L)
  expect_equal(cigar_ref_width("20M5I20M"), 40L)
  expect_equal(cigar_ref_width("10M100N10M"), 120L)
  expect_equal(cigar_ref_width(c("50M", "10S40M")), c(50L, 40L))
  expect_true(is.na(cigar_ref_width("*")))
})

# reads: one plain, one clipped on each side, one with a deletion in the
# CIGAR, and a locus at 700-710 with both intact and clipped reads so the
# genotype call has something to separate.
mini_bam <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dir <- fs::path(tempdir(), "wca_mini_bam")
    fs::dir_create(dir)
    sam <- fs::path(dir, "mini.sam")
    a50 <- strrep("A", 50)
    polya <- strrep("A", 20)
    hexa <- strrep("CCCTCT", 4)   # 24 bases
    body <- c(
      # qname flag rname pos mapq cigar rnext pnext tlen seq qual
      str_c("r1\t99\tchrT\t100\t60\t50M\t=\t300\t250\t", a50, "\t", strrep("I", 50)),
      str_c("r2\t99\tchrT\t200\t60\t10S40M\t=\t400\t240\t", strrep("T", 10), strrep("G", 40), "\t", strrep("I", 50)),
      str_c("r3\t99\tchrT\t300\t60\t40M10S\t=\t500\t240\t", strrep("G", 40), strrep("C", 10), "\t", strrep("I", 50)),
      str_c("r4\t99\tchrT\t400\t60\t20M5D20M\t=\t600\t240\t", strrep("G", 40), "\t", strrep("I", 40)),
      # at the locus 700-710: two intact spanning reads
      str_c("r5\t99\tchrT\t680\t60\t60M\t=\t900\t280\t", strrep("G", 60), "\t", strrep("I", 60)),
      str_c("r6\t99\tchrT\t685\t60\t60M\t=\t905\t280\t", strrep("G", 60), "\t", strrep("I", 60)),
      # one clipped left at 700 carrying poly-A, one clipped right at 710
      # carrying the hexamer
      str_c("r7\t99\tchrT\t700\t60\t20S40M\t=\t950\t290\t", polya, strrep("G", 40), "\t", strrep("I", 60)),
      str_c("r8\t99\tchrT\t671\t60\t40M24S\t=\t400\t-290\t", strrep("G", 40), hexa, "\t", strrep("I", 64))
    )
    write_lines(c("@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chrT\tLN:2000", body), sam)
    bam <- Rsamtools::asBam(sam, fs::path(dir, "mini"), overwrite = TRUE,
                            indexDestination = TRUE)
    cache <<- as.character(bam)
    cache
  }
})

test_that("bam_reads reports the aligned span, not the read span", {
  skip_without_bam()
  r <- bam_reads(mini_bam(), "chrT", 1, 2000)
  expect_true(all(c("read_id", "pos", "ref_end", "mapq", "cigar", "seq") %in% names(r)))
  expect_equal(nrow(r), 8L)

  plain <- r |> filter(cigar == "50M")
  expect_equal(plain$pos, 100L)
  expect_equal(plain$ref_end, 149L)

  # soft clips consume no reference, so a 10S40M read spans 40 bases
  clipped <- r |> filter(cigar == "10S40M")
  expect_equal(clipped$pos, 200L)
  expect_equal(clipped$ref_end, 239L)

  # a deletion does consume reference
  deleted <- r |> filter(cigar == "20M5D20M")
  expect_equal(deleted$ref_end, 444L)
})

test_that("bam_reads returns an empty typed tibble for an empty window", {
  skip_without_bam()
  r <- bam_reads(mini_bam(), "chrT", 1800, 1900)
  expect_equal(nrow(r), 0L)
  expect_true(is.integer(r$pos))
  expect_true(is.character(r$cigar))
})

test_that("bam_soft_clips places each clip at the reference base it abuts", {
  skip_without_bam()
  cl <- bam_reads(mini_bam(), "chrT", 1, 2000) |> bam_soft_clips(min_len = 5)

  left <- cl |> filter(side == "left", clip_len == 10L)
  expect_equal(left$clip_pos, 200L)          # first aligned base
  expect_equal(left$clip_seq, strrep("T", 10))

  right <- cl |> filter(side == "right", clip_len == 10L)
  expect_equal(right$clip_pos, 339L)         # last aligned base
  expect_equal(right$clip_seq, strrep("C", 10))
  expect_equal(right$gc_frac, 1)

  # an unclipped read contributes no rows
  expect_false(any(cl$clip_len == 0))
  expect_equal(sum(cl$at_frac == 1), 2L)     # the poly-T and the poly-A clip
})

test_that("bam_soft_clips honours min_len", {
  skip_without_bam()
  r <- bam_reads(mini_bam(), "chrT", 1, 2000)
  expect_equal(nrow(bam_soft_clips(r, min_len = 15L)), 2L)   # the 20S and the 24S
  expect_equal(nrow(bam_soft_clips(r, min_len = 100L)), 0L)
})

test_that("bam_base_depth walks the CIGAR, so a deletion reads as uncovered", {
  skip_without_bam()
  d <- bam_base_depth(mini_bam(), "chrT", 395, 450)
  expect_equal(nrow(d), 56L)
  expect_equal(d$depth[d$pos == 399], 0L)    # before the read
  expect_equal(d$depth[d$pos == 400], 1L)    # first aligned base
  expect_equal(d$depth[d$pos == 419], 1L)    # last base of the first block
  expect_equal(d$depth[d$pos == 420], 0L)    # inside the 5 base deletion
  expect_equal(d$depth[d$pos == 424], 0L)
  expect_equal(d$depth[d$pos == 425], 1L)    # second block resumes
  expect_equal(d$depth[d$pos == 444], 1L)
  expect_equal(d$depth[d$pos == 445], 0L)

  # soft-clipped bases are not coverage either
  d2 <- bam_base_depth(mini_bam(), "chrT", 190, 205)
  expect_equal(d2$depth[d2$pos == 195], 0L)
  expect_equal(d2$depth[d2$pos == 200], 1L)
})

test_that("bam_binned_counts is relative to the median bin", {
  skip_without_bam()
  b <- bam_binned_counts(mini_bam(), "chrT", 100, 1000, n_bins = 9L, bin_width = 100L)
  expect_equal(nrow(b), 9L)
  expect_equal(sum(b$n > 0) > 0, TRUE)
  expect_equal(b$bin_end, b$bin_start + 99L)
  expect_equal(b$rel, b$n / stats::median(b$n))
})

test_that("clip_genotype separates intact reads from element clips", {
  skip_without_bam()
  classes <- list(
    hexamer = \(s) str_detect(s, "CCCTCT"),
    polyA = \(s) !str_detect(s, "CCCTCT") & str_count(s, "[AT]") / nchar(s) > 0.9)
  g <- clip_genotype(mini_bam(), "chrT", 700, 710, classes = classes,
                     window = 200L, min_clip = 15L, flank = 5L)
  expect_equal(g$n_intact, 2L)       # r5 and r6 cross 695-715 with no clip
  expect_equal(g$n_clipped, 2L)
  expect_equal(g$hexamer, 1L)
  expect_equal(g$polyA, 1L)
  expect_equal(g$vaf_clipped, 0.5)
})

test_that("clip_genotype at a locus with no clips calls the reference allele", {
  skip_without_bam()
  g <- clip_genotype(mini_bam(), "chrT", 120, 130,
                     classes = list(anything = \(s) rep(TRUE, length(s))),
                     window = 50L, min_clip = 15L, flank = 5L)
  expect_equal(g$n_clipped, 0L)
  expect_equal(g$anything, 0L)
  expect_equal(g$n_intact, 1L)
  expect_equal(g$vaf_clipped, 0)
})

test_that("ref_at_profile finds a pure A/T tract", {
  testthat::skip_if_not_installed("Rsamtools")
  dir <- fs::path(tempdir(), "wca_mini_fa")
  fs::dir_create(dir)
  fa <- fs::path(dir, "mini.fa")
  # 100 bases of mixed sequence, then 100 of pure A/T, then mixed again
  mixed <- strrep("ACGT", 25)
  tract <- strrep("TAAA", 25)
  write_lines(c(">chrT", str_c(mixed, tract, mixed)), fa)
  Rsamtools::indexFa(fa)

  p <- ref_at_profile(fa, "chrT", 1, 300, window = 20L, step = 10L)
  expect_true(all(c("chrom", "pos", "at_frac") %in% names(p)))
  # inside the tract every window is all A and T
  inside <- p |> filter(pos > 120, pos < 190)
  expect_true(all(inside$at_frac == 1))
  # inside the mixed sequence exactly half
  outside <- p |> filter(pos > 20, pos < 80)
  expect_true(all(outside$at_frac == 0.5))
})

test_that("read_pair_orientation reads FR, RF, FF and RR off the flag", {
  # flag bits: 16 = this read reverse, 32 = mate reverse
  reads <- tibble(
    flag = c(99L, 16L + 1L, 1L, 16L + 32L + 1L, 32L + 1L),
    pos = c(100L, 100L, 100L, 100L, 500L),
    mate_pos = c(500L, 500L, 500L, 500L, 100L))
  o <- read_pair_orientation(reads)
  # 99 = paired, proper, mate reverse: left mate forward, right reverse
  expect_equal(o$orientation[1], "FR")
  # self reverse, mate forward, self on the left
  expect_equal(o$orientation[2], "RF")
  # neither reverse
  expect_equal(o$orientation[3], "FF")
  # both reverse
  expect_equal(o$orientation[4], "RR")
  # self on the right and reverse-flagged mate: orientation is read from
  # the leftmost mate, so this is FR again, not RF
  expect_equal(o$orientation[5], "RF")
  expect_equal(o$sv_signature[1], "deletion or normal")
  expect_equal(o$sv_signature[3], "inversion, head to head")
  expect_equal(o$sv_signature[4], "inversion, tail to tail")
})

test_that("read_pair_orientation leaves an unmapped mate as NA", {
  o <- read_pair_orientation(tibble(flag = 1L, pos = 100L, mate_pos = NA_integer_))
  expect_true(is.na(o$orientation))
  expect_true(is.na(o$sv_signature))
})

test_that("clip_partner_match tests both orientations of the partner locus", {
  partner <- "GGGGCCCCATCGATCGATCGGATCCAAGCTTGGGGCCCCTTTTAAAA"
  # a clip that is literally partner sequence: a deletion junction
  fwd_clip <- "ATCGATCGATCGGATCC"
  # a clip that is the reverse complement: an inversion junction
  rev_clip <- revcomp(fwd_clip)
  # a clip from neither: an insertion of sequence not in the reference
  alien_clip <- strrep("CCCTCT", 6)

  clips <- tibble(clip_seq = c(fwd_clip, rev_clip, alien_clip),
                  clip_len = nchar(c(fwd_clip, rev_clip, alien_clip)))
  m <- clip_partner_match(clips, partner)

  expect_equal(m$match_forward[1], nchar(fwd_clip))
  expect_true(m$match_revcomp[1] < nchar(fwd_clip))
  expect_equal(m$match_frac[1], 1)

  expect_equal(m$match_revcomp[2], nchar(rev_clip))
  expect_true(m$match_forward[2] < nchar(rev_clip))
  expect_equal(m$match_frac[2], 1)

  # the alien clip matches neither orientation beyond chance
  expect_true(m$match_frac[3] < 0.3)
})

test_that("clip_partner_match refuses an empty partner sequence", {
  clips <- tibble(clip_seq = "ACGT", clip_len = 4L)
  expect_error(clip_partner_match(clips, NA_character_), "partner_seq is empty")
  expect_error(clip_partner_match(clips, ""), "partner_seq is empty")
})

test_that("revcomp round trips and complements", {
  expect_equal(revcomp("ACGT"), "ACGT")
  expect_equal(revcomp("AAAA"), "TTTT")
  expect_equal(revcomp("CCCTCT"), "AGAGGG")
  expect_equal(revcomp(revcomp("ACGTTGCA")), "ACGTTGCA")
  expect_equal(revcomp(c("AC", "GT")), c("GT", "AC"))
})
