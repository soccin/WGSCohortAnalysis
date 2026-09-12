mini_counts <- function() {
  dir <- fs::path(tempdir(), "wca_mini_fc")
  fs::dir_create(dir)
  f <- fs::path(dir, "mini.gene.featureCounts.txt")
  write_lines(c(
    "# Program:featureCounts v2.0.1; Command:whatever",
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\tsome.sample.bam",
    "ENSG00000000001\t1\t1\t1000\t+\t1000\t500",
    "ENSG00000000002\t1\t2000\t4000\t-\t2000\t0",
    "ENSG00000000003\t2\t1\t500\t+\t500\t500"), f)
  f
}

test_that("read_featurecounts_genes computes CPM and FPKM from the whole library", {
  d <- read_featurecounts_genes(mini_counts(),
    c(GENE_A = "ENSG00000000001", GENE_B = "ENSG00000000002"))
  expect_equal(d$gene, c("GENE_A", "GENE_B"))
  expect_equal(d$count, c(500, 0))
  # library size is every gene in the file, not just the two asked for
  expect_equal(unique(d$library_size), 1000)
  expect_equal(d$CPM, c(500000, 0))
  # 500 counts / (1000/1e6) per million / (1000/1e3) kb
  expect_equal(d$FPKM[1], 500 / (1000 / 1e6) / 1)
})

test_that("a gene missing from the file comes back at zero, not as NA", {
  d <- read_featurecounts_genes(mini_counts(),
    c(GENE_A = "ENSG00000000001", NOT_THERE = "ENSG09999999999"))
  expect_equal(d$count, c(500, 0))
  expect_equal(d$CPM[2], 0)
  expect_true(is.na(d$FPKM[2]))
})

test_that("read_featurecounts_genes refuses an empty gene list or a missing file", {
  expect_error(read_featurecounts_genes(mini_counts(), character()), "gene_ids is empty")
  expect_error(read_featurecounts_genes(fs::path(tempdir(), "nope.txt"),
                                        c(A = "ENSG1")), "not found")
})

test_that("featurecounts_beside finds the counts next to a fusion call", {
  root <- fs::path(tempdir(), "wca_mini_forte", "SAMPLE_A")
  fs::dir_create(fs::path(root, "metafusion"))
  fs::dir_create(fs::path(root, "featurecounts"))
  cff <- fs::path(root, "metafusion", "SAMPLE_A.final.cff")
  counts <- fs::path(root, "featurecounts", "SAMPLE_A.gene.featureCounts.txt")
  write_lines("", cff)
  write_lines("", counts)
  expect_equal(featurecounts_beside(cff), as.character(counts))

  lonely <- fs::path(tempdir(), "wca_mini_forte", "SAMPLE_B", "metafusion")
  fs::dir_create(lonely)
  orphan <- fs::path(lonely, "SAMPLE_B.final.cff")
  write_lines("", orphan)
  expect_true(is.na(featurecounts_beside(orphan)))
})
