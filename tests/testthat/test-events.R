events_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      coh <- fixture_cohort()
      maf <- read_cohort_maf(coh)
      fac <- read_cohort_facets(coh, drop_diploid = FALSE)
      sv_raw <- read_cohort_sv(coh)
      qc <- facets_qc_gate(fac$qc)
      snv <- snv_events(maf)
      cnv <- cnv_events(fac$gene, qc)
      arm <- cnv_arm_events(fac$arm, qc)
      sv <- sv_events(sv_raw)
      svg <- sv_gene_events(sv)
      cache <<- list(coh = coh, maf = maf, fac = fac, qc = qc, snv = snv, cnv = cnv, arm = arm,
        sv = sv, svg = svg, events = bind_events(snv, cnv, svg))
    }
    cache
  }
})

test_that("snv_events keeps coding rows and TERT promoter, flags LoF", {
  e <- events_fixture()
  snv <- e$snv
  expect_true(all(!is.na(snv$Alteration) | snv$VClass == "5'Flank"))
  expect_true(any(snv$Gene == "TERT" & snv$VClass == "5'Flank"))
  expect_true(any(snv$is_lof))
  expect_true(all(snv$is_lof == (snv$VClass %in% maf_lof_classes())))
  expect_equal(names(snv)[1:3], c("Gene", "Sample", "NID"))
  strict <- snv_events(e$maf, min_t_depth = 1e6)
  expect_equal(nrow(strict), 0)
})

test_that("cnv_events drops diploid, X and excluded samples, and calls Dir", {
  e <- events_fixture()
  cnv <- e$cnv
  expect_false(any(cnv$cn_state == "DIPLOID"))
  expect_true(all(cnv$chrom %in% autosomes()))
  expect_true(all(cnv$Dir %in% c("gain", "loss", "cnloh", NA)))
  with_x <- cnv_events(e$fac$gene, e$qc, autosomes_only = FALSE)
  expect_true("X" %in% with_x$chrom)
  qc_all_excluded <- mutate(e$qc, EXCLUDE = TRUE)
  expect_equal(nrow(cnv_events(e$fac$gene, qc_all_excluded)), 0)
  tcn_rule <- cnv_events(e$fac$gene, e$qc, diploid_rule = "tcn")
  expect_false(any(tcn_rule$tcn == 2))
})

test_that("cnv_arm_events is autosomal with a Dir", {
  arm <- events_fixture()$arm
  expect_true(all(arm$chrom %in% autosomes()))
  expect_true("Dir" %in% names(arm))
})

test_that("every CNV summarizer and plot drops X by default even on unfiltered input", {
  e <- events_fixture()
  with_x <- cnv_events(e$fac$gene, e$qc, autosomes_only = FALSE)
  arm_x <- cnv_arm_events(e$fac$arm, e$qc, autosomes_only = FALSE)
  expect_true("X" %in% with_x$chrom)
  expect_true("X" %in% arm_x$chrom)
  expect_true(all(filter_autosomes(with_x)$chrom %in% autosomes()))
  expect_equal(nrow(filter_autosomes(tibble(chrom = c("23", "chrX", "Y", "7", "chr1")))), 2)

  bands <- cytoband_recurrence(with_x, 3)
  expect_false(any(bands$chrom %in% c("X", "Y")))
  expect_equal(bands, cytoband_recurrence(e$cnv, 3))
  bands_x <- cytoband_recurrence(with_x, 3, autosomes_only = FALSE)
  expect_true("X" %in% bands_x$chrom)

  arms <- arm_recurrence(arm_x, 3)
  expect_false(any(arms$chrom %in% c("X", "Y")))
  expect_equal(arms, arm_recurrence(e$arm, 3))
  expect_true("X" %in% arm_recurrence(arm_x, 3, autosomes_only = FALSE)$chrom)

  coh <- fixture_cohort()
  s <- sample_summary(coh, cnv_gene = with_x)
  expect_equal(s$n_cnv_genes, sample_summary(coh, cnv_gene = e$cnv)$n_cnv_genes)
  expect_true(any(sample_summary(coh, cnv_gene = with_x, autosomes_only = FALSE)$n_cnv_genes > s$n_cnv_genes))

  expect_equal(plot_cytoband_gainloss(bands_x)$data, plot_cytoband_gainloss(bands)$data)
  arm_levels <- levels(plot_arm_heatmap(arm_x)$data$arm)
  expect_false(any(str_detect(arm_levels, "^[XY]")))
  expect_equal(arm_levels, levels(plot_arm_heatmap(e$arm)$data$arm))
  expect_true(any(str_detect(levels(plot_arm_heatmap(arm_x, autosomes_only = FALSE)$data$arm), "^X")))
})

test_that("sv_events annotates bands, pair keys and fusion classes", {
  sv <- events_fixture()$sv
  expect_true(all(!is.na(sv$bandA)))
  expect_true(all(sv$PairKey == gene_pair_key(sv$gene1, sv$gene2)))
  expect_true("in-frame" %in% sv$fusion_class)
  expect_true("intragenic" %in% sv$fusion_class)
  expect_true(all(sv$fusion_class[sv$fusion == "-"] == "none"))
  strict <- sv_events(events_fixture()$sv |> select(-PairKey, -fusion_class, -is_fusion, -bandA, -bandB), min_callers_pass = 3)
  expect_true(all(strict$NumCallersPass >= 3))
})

test_that("gene_pair_key is order independent", {
  expect_equal(gene_pair_key("B", "A"), "A::B")
  expect_equal(gene_pair_key(c("A", "Z"), c("B", "C")), c("A::B", "C::Z"))
})

test_that("sv_gene_events splits pairs and dedups per Gene-UUID", {
  e <- events_fixture()
  svg <- e$svg
  expect_false(any(duplicated(svg[, c("Gene", "UUID")])))
  self <- e$sv |> filter(gene1 == gene2)
  expect_gt(nrow(self), 0)
  expect_equal(sum(svg$UUID %in% self$UUID), n_distinct(self$UUID))
  inter <- e$sv |> filter(gene1 != gene2)
  expect_equal(sum(svg$UUID %in% inter$UUID), 2 * n_distinct(inter$UUID))
})

test_that("bind_events produces the long table", {
  ev <- events_fixture()$events
  expect_equal(names(ev), c("Sample", "Gene", "EventType", "Detail"))
  expect_true(all(ev$EventType %in% c("SNV", "SNV-LoF", "CNV-gain", "CNV-loss", "CNV-cnloh", "SV", "Fusion")))
  expect_true(all(c("SNV", "SV", "CNV-loss") %in% ev$EventType))
})

test_that("recurrence_table honors the denominator", {
  e <- events_fixture()
  r <- recurrence_table(e$snv, 10, "SNV")
  expect_equal(attr(r, "denominator"), 10)
  expect_equal(r$PCT, r$n / 10)
  expect_equal(names(r), c("Gene", "n", "PCT", "Event"))
  expect_true(all(r$n <= 3))
  expect_error(recurrence_table(e$snv, event = "SNV"), "denominator")
})

test_that("event_summary has the expected column set and PCT", {
  e <- events_fixture()
  tabs <- list(snv = recurrence_table(e$snv, 3, "SNV"), cnv = recurrence_table(e$cnv, 3, "CNV"),
    sv = recurrence_table(e$svg, 3, "SV"))
  s <- event_summary(tabs, 3)
  expect_equal(names(s), c("Gene", "Nt", "PCT", "CNV_n", "SNV_n", "SV_n", "CNV_PCT", "SNV_PCT", "SV_PCT"))
  expect_equal(s$Nt, s$CNV_n + s$SNV_n + s$SV_n)
  expect_equal(s$PCT, s$Nt / 3)
  expect_true(all(diff(s$Nt) <= 0))
  expect_equal(gene_set_table(s, c("TP53", "GATA3"))$Gene, sort(intersect(c("TP53", "GATA3"), s$Gene)))
})

test_that("cytoband, arm, pair, partner and fusion tables", {
  e <- events_fixture()
  bands <- cytoband_recurrence(e$cnv, 3)
  expect_equal(names(bands)[1:7], c("Band", "chrom", "N_total", "N_gain", "PCT_gain", "N_loss", "PCT_loss"))
  expect_true(all(bands$N_total >= bands$N_gain + bands$N_loss))
  expect_true(all(bands$N_gain <= 3))
  arms <- arm_recurrence(e$arm, 3)
  expect_equal(nrow(arms), nrow(filter_autosomes(chrom_arms())))
  expect_equal(nrow(arm_recurrence(e$arm, 3, autosomes_only = FALSE)), nrow(chrom_arms()))
  expect_true(all(arms$PCT_loss <= 1))
  pairs <- sv_pair_recurrence(e$sv, 3)
  expect_true(all(pairs$n <= 3))
  expect_equal(pairs$PairKey, gene_pair_key(pairs$gene1, pairs$gene2))
  partners <- sv_partner_table(e$svg)
  expect_equal(names(partners), c("Gene", "N", "GenePartners", "Partners"))
  fus <- fusion_table(e$sv, 3)
  expect_true(nrow(fus) > 0)
  expect_true(all(fus$fusion_class %in% c("in-frame", "out-of-frame", "protein-fusion", "transcript")))
})

test_that("group_recurrence counts per group and ignores ungrouped samples", {
  e <- events_fixture()
  groups <- list(CellLine = c("CL01", "CL02"), Patient = "CL03")
  g <- group_recurrence(e$svg, groups, detail_cols = "sv_class")
  expect_equal(names(g), c("Gene", "n", "n_CellLine", "PCT_CellLine", "Samples_CellLine",
    "n_Patient", "PCT_Patient", "Samples_Patient", "sv_class"))
  expect_equal(attr(g, "denominator"), c(CellLine = 2L, Patient = 1L))
  expect_equal(g$n, g$n_CellLine + g$n_Patient)
  expect_equal(g$PCT_CellLine, g$n_CellLine / 2)
  expect_true(all(diff(g$n) <= 0))
  expect_true(all(str_detect(g$Samples_Patient, "^(CL03)?$")))
  both <- g |> filter(n_CellLine == 2, n_Patient >= 1)
  expect_true(all(both$Samples_CellLine == "CL01;CL02"))
  expect_setequal(both$Gene, e$svg |> distinct(Gene, Sample) |> count(Gene) |> filter(n == 3) |> pull(Gene))
  # ungrouped samples drop out of n
  one <- group_recurrence(e$svg, list(A = "CL01"))
  expect_true(all(one$n == 1))
  expect_equal(nrow(one), n_distinct(filter(e$svg, Sample == "CL01")$Gene))
  # pair-level keys work too
  p <- group_recurrence(e$sv, groups, key_col = "PairKey", detail_cols = c("sv_class", "fusion_class"))
  expect_equal(names(p)[1], "PairKey")
  expect_true(all(p$n <= 3))
  expect_error(group_recurrence(e$svg, list(c("CL01"))), "named")
  expect_error(group_recurrence(e$svg, list(A = "CL01", B = "CL01")), "more than one group")
  expect_error(group_recurrence(e$svg, list(`Cell Line` = "CL01")), "alphanumeric")
  d <- group_recurrence_dictionary(groups)
  expect_equal(d$Column, c("n_CellLine", "PCT_CellLine", "Samples_CellLine", "n_Patient", "PCT_Patient", "Samples_Patient"))
})

test_that("sample_summary has one row per sample and blanks failed QC", {
  e <- events_fixture()
  qcd <- read_cohort_qc(e$coh)
  s <- sample_summary(e$coh, snv = e$snv, sv = e$sv, cnv_gene = e$cnv, facets_qc = e$qc,
    sample_data = qcd$sample_data, qc_status = qcd$qc_status)
  expect_equal(nrow(s), 3)
  expect_true(all(c("n_snv", "n_sv", "n_sv_TRA", "TMB", "purity", "QC_Status", "n_cnv_genes") %in% names(s)))
  expect_equal(s$QC_Status, rep("pass", 3))
  qc_fail <- mutate(e$qc, facets_qc = FALSE)
  s2 <- sample_summary(e$coh, facets_qc = qc_fail)
  expect_true(all(is.na(s2$purity)))
})

test_that("sv_to_bed yields two rows per SV", {
  sv <- events_fixture()$sv
  bed <- sv_to_bed(sv)
  expect_equal(nrow(bed), 2 * nrow(sv))
  expect_equal(names(bed)[1:6], c("chrom", "start", "end", "name", "score", "strand"))
  f <- withr::local_tempfile(fileext = ".bed")
  write_sv_bed(sv, f)
  expect_true(str_starts(readr::read_lines(f, n_max = 1), "track"))
})
