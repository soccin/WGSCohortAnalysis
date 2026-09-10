test_that("normalize_chrom handles prefixes and FACETS codes", {
  expect_equal(normalize_chrom(c("chr23", "24", "1", "chrX", "M", 25, "chr7")),
    c("X", "Y", "1", "X", "MT", "MT", "7"))
  expect_equal(levels(chrom_factor(c("2", "X", "chr1"))), chrom_levels("hg19"))
  expect_equal(as.character(sort(chrom_factor(c("X", "10", "2")))), c("2", "10", "X"))
})

test_that("genome_info derives lengths, centromeres and arms from cytoBand", {
  g <- genome_info("hg19")
  expect_equal(g$chrom_levels, c(as.character(1:22), "X", "Y"))
  expect_equal(unname(g$chrom_len[["1"]]), 249250621)
  expect_equal(unname(g$chrom_len[["X"]]), 155270560)
  cen1 <- g$centromeres |> filter(chrom == "1")
  expect_equal(cen1$cen_start, 121500000)
  expect_equal(cen1$cen_end, 128900000)
  arms <- chrom_arms("hg19")
  expect_false(any(arms$arm %in% c("13p", "14p", "15p", "21p", "22p")))
  expect_true(all(c("13q", "17p", "17q", "Xp", "Xq") %in% arms$arm))
  expect_equal(arms |> filter(arm == "17p") |> pull(end), 24000000)
  expect_equal(autosomes(), as.character(1:22))
  expect_true("GRCh37" %in% g$build_aliases)
})

test_that("annotate_cytoband joins a known gene position", {
  btk <- tibble(chrom = "X", gene_start = 100604435)
  expect_equal(annotate_cytoband(btk)$Band, "Xq22.1")
  tp53 <- tibble(chrom = "chr17", pos = 7571720)
  expect_equal(annotate_cytoband(tp53, "chrom", "pos", name = "band")$band, "17p13.1")
  expect_equal(annotate_arm(tp53, "chrom", "pos")$arm, "17p")
  off <- tibble(chrom = "1", gene_start = 9e9)
  expect_true(is.na(annotate_cytoband(off)$Band))
})

test_that("unknown genomes are rejected", {
  expect_error(genome_info("hg38"), "Unknown genome")
})
