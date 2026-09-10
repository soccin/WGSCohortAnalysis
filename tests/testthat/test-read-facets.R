test_that("tempo_pair_files finds every expected file", {
  pf <- tempo_pair_files(fixture_pair("CL01"))
  expect_equal(pf$Sample, "CL01")
  expect_equal(pf$NID, "CL01N")
  for (col in c("maf", "bedpe", "facets_gene", "facets_arm", "facets_qc", "facets_out",
    "facets_cncf", "facets_seg", "facets_fit_qc", "sample_data", "qc_status")) {
    expect_true(fs::file_exists(pf[[col]]), info = col)
  }
  expect_true(is.na(pf$maf_unfiltered))
  expect_equal(tempo_pair_files(pf$maf)$pair_dir, pf$pair_dir)
  cf <- tempo_cohort_files(pf$maf)
  expect_true(fs::file_exists(cf$sample_data))
  expect_true(fs::file_exists(cf$alignment_qc))
})

test_that("read_facets_gene maps chrom 23 to X and splits the pair name", {
  pf <- tempo_pair_files(fixture_pair("CL01"))
  g <- read_facets_gene(pf$facets_gene)
  expect_true("X" %in% g$chrom)
  expect_false("23" %in% g$chrom)
  expect_equal(unique(g$Sample), "CL01")
  expect_equal(unique(g$NID), "CL01N")
  g2 <- read_facets_gene(pf$facets_gene, drop_diploid = TRUE)
  expect_false(any(g2$cn_state == "DIPLOID"))
})

test_that("read_facets_arm, qc, cncf and seg parse", {
  pf <- tempo_pair_files(fixture_pair("CL02"))
  arm <- read_facets_arm(pf$facets_arm)
  expect_true(all(c("Sample", "arm", "chrom", "frac_of_arm", "cn_state") %in% names(arm)))
  expect_true(is.numeric(arm$frac_of_arm))
  qc <- read_facets_qc(pf$facets_qc)
  expect_equal(nrow(qc), 1)
  expect_true(is.logical(qc$facets_qc))
  expect_true(is.numeric(qc$dipLogR))
  cncf <- read_facets_cncf(pf$facets_cncf)
  expect_equal(unique(cncf$Sample), "CL02")
  seg <- read_facets_seg(pf$facets_seg)
  expect_true(all(c("chrom", "loc.start", "loc.end", "seg.mean") %in% names(seg)))
})

test_that("cn_call presets", {
  st <- c("HETLOSS", "HOMDEL", "GAIN", "AMP", "CNLOH", "DIPLOID", "LOSS BEFORE", "AMP (LOH)", "CNLOH & GAIN")
  expect_equal(cn_call(st, preset = "facets"),
    c("loss", "loss", "gain", "gain", "cnloh", NA, "loss", "gain", "gain"))
  expect_equal(cn_call(st, preset = "strict"),
    c("loss", "loss", "gain", NA, NA, NA, NA, NA, NA))
  expect_equal(cn_call(NULL, tcn = c(1, 3, 2, 2), lcn = c(0, 1, 0, 1), preset = "tcn"),
    c("loss", "gain", "cnloh", NA))
  expect_equal(cn_call(NULL, tcn = c(1, 1, 3), lcn = c(0, 0, 1), filter = c("PASS", "FAIL", "RESCUE"),
    genes_on_seg = c(3, 3, 50), preset = "focal"), c("loss", NA, NA))
})

test_that("facets_qc_gate excludes on dipLogR and optionally on facets_qc", {
  qc <- tibble(Sample = c("a", "b", "c"), facets_qc = c(TRUE, FALSE, TRUE), dipLogR = c(0.1, 0.2, -2))
  g <- facets_qc_gate(qc)
  expect_equal(g$EXCLUDE, c(FALSE, FALSE, TRUE))
  expect_equal(cnv_samples(g), c("a", "b"))
  g2 <- facets_qc_gate(qc, require_facets_qc = TRUE)
  expect_equal(g2$EXCLUDE, c(FALSE, TRUE, TRUE))
  expect_equal(g2$exclude_reason[2], "facets_qc FALSE")
})

test_that("read_cohort_facets returns gene, arm and qc keyed by cohort ids", {
  fac <- read_cohort_facets(fixture_cohort())
  expect_setequal(names(fac), c("gene", "arm", "qc"))
  expect_equal(nrow(fac$qc), 3)
  expect_setequal(unique(fac$gene$Sample), c("CL01", "CL02", "CL03"))
  expect_false("Sample_file" %in% names(fac$gene))
})

test_that("QC and sample_data readers", {
  pf <- tempo_pair_files(fixture_pair("CL03"))
  qs <- read_qc_status(pf$qc_status)
  expect_equal(names(qs), c("Sample", "Status", "Reason"))
  expect_equal(qs$Sample, c("CL03", "CL03N"))
  sd <- read_sample_data(pf$sample_data)
  expect_equal(sd$Sample, "CL03")
  expect_true(is.numeric(sd$TMB))
  expect_true(is.logical(sd$WGD_status))
  expect_false(any(str_starts(names(sd), "SBS")))
  aq <- read_alignment_qc(tempo_cohort_files(pf$maf)$alignment_qc)
  expect_equal(nrow(aq), 6)
  qc <- read_cohort_qc(fixture_cohort())
  expect_equal(sum(qc$qc_status$is_tumor), 3)
  expect_equal(nrow(qc$sample_data), 3)
})
