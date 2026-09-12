test_that("longest_common_substring finds the match and its position", {
  r <- longest_common_substring("AACCGGTTAA", "ZZCCGGTTZZ")
  expect_equal(r$len, 6L)
  expect_equal(r$seq, "CCGGTT")
  expect_equal(r$a_start, 3L)
  expect_equal(r$b_start, 3L)

  expect_equal(longest_common_substring("ACGT", "ZZZZ")$len, 0L)
  expect_equal(longest_common_substring("", "ACGT")$len, 0L)
  expect_equal(longest_common_substring("ACGT", "ACGT")$len, 4L)
})

test_that("junction_anchor_support separates a unique anchor from an A/T tract", {
  # A consensus built from a poly-A tract joined to unique sequence, the
  # shape a mobile element insertion produces when the caller reports it
  # as a deletion between two loci.
  tract <- strrep("A", 60)
  unique_seq <- "GCCTAGGCATCGGATCCTAGCATCGGATCCTAGCATCGGATCCTAGCATCGG"
  consensus <- str_c(tract, unique_seq)

  loci <- tibble(
    breakend = c("A", "B"),
    chrom = c("6", "6"),
    pos = c(1000L, 2000L),
    seq = c(str_c("TTTT", tract, "TTTT"), str_c("TTTT", unique_seq, "TTTT"))
  )
  out <- junction_anchor_support(consensus, loci)

  expect_equal(nrow(out), 2)
  expect_equal(out$match_len[out$breakend == "A"], 60L)
  expect_equal(out$gc_in_match[out$breakend == "A"], 0L)
  expect_false(out$informative[out$breakend == "A"])
  expect_equal(out$at_fraction[out$breakend == "A"], 1)

  expect_equal(out$match_len[out$breakend == "B"], str_length(unique_seq))
  expect_gt(out$gc_in_match[out$breakend == "B"], 5)
  expect_true(out$informative[out$breakend == "B"])

  expect_error(junction_anchor_support(consensus, loci |> select(-seq)), "loci\\$seq or fasta")
  expect_error(junction_anchor_support(consensus, tibble(breakend = "A")), "missing columns")
})

test_that("plot_junction_anatomy draws one row per breakend", {
  loci <- tibble(breakend = c("A", "B"), chrom = "6", pos = c(1000L, 2000L),
                 seq = c(strrep("A", 40), "GCCTAGGCATCGGATCCTAGCATCGG"))
  out <- junction_anchor_support(str_c(strrep("A", 40), "GCCTAGGCATCGGATCCTAGCATCGG"), loci) |>
    mutate(Sample = breakend)
  p <- plot_junction_anatomy(out, consensus_len = 66L, title = "test")
  expect_s3_class(p, "ggplot")

  # A consensus_len column in the data must not mask the argument: the
  # background segment silently lost its end point when it did.
  with_col <- out |> mutate(consensus_len = 66L)
  lens <- set_names(c(66L, 66L), with_col$Sample)
  built <- ggplot2::ggplot_build(plot_junction_anatomy(with_col, consensus_len = lens))
  expect_false(any(is.na(built$data[[1]]$xend)))

  # A single length applies to every row; a named vector has to cover them.
  expect_s3_class(plot_junction_anatomy(with_col, consensus_len = 66L), "ggplot")
  expect_error(plot_junction_anatomy(with_col, consensus_len = c(nope = 66L, other = 66L)),
               "no consensus length")
})
