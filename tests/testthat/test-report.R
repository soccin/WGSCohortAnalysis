test_that("write_report_xlsx writes sheets, PCT formats and a DataDictionary", {
  t1 <- tibble(Gene = c("A", "B"), n = c(2, 1), PCT = c(0.5, 0.25), Event = "SNV")
  attr(t1, "denominator") <- 4
  t2 <- tibble(Sample = "s1", weird_col = 1)
  f <- withr::local_tempfile(fileext = ".xlsx")
  wb <- write_report_xlsx(list(SNVGenes = t1, Other = t2), f, notes = c(cohort = "4 samples"))
  expect_true(fs::file_exists(f))
  expect_equal(readxl::excel_sheets(f), c("SNVGenes", "Other", "DataDictionary"))
  back <- readxl::read_excel(f, "SNVGenes")
  expect_equal(back$PCT, c(0.5, 0.25))
  dict <- readxl::read_excel(f, "DataDictionary")
  expect_true(any(dict$Sheet == "(notes)" & dict$Column == "cohort"))
  expect_true(any(dict$Sheet == "SNVGenes" & dict$Column == "(denominator)" & str_detect(dict$Description, "4 samples")))
  expect_true(any(dict$Sheet == "Other" & dict$Column == "weird_col"))
  styles <- openxlsx::getStyles(wb)
  expect_true(any(purrr::map_lgl(styles, \(s) identical(s$numFmt$formatCode, "0.00%"))))
  expect_error(write_report_xlsx(list(t1), f), "named")
})

test_that("wca_dictionary is unique per sheet and column", {
  d <- wca_dictionary()
  expect_false(any(duplicated(d[, c("Sheet", "Column")])))
})

test_that("cached_read keys on content, reader version and arguments", {
  dir <- withr::local_tempdir()
  f <- fs::path(dir, "x.txt")
  readr::write_lines(c("a", "b"), f)
  reader <- function(file, n = 1) readr::read_lines(file, n_max = n)
  expect_equal(cached_read(f, reader, cache_dir = dir), "a")
  expect_equal(cached_read(f, reader, n = 2, cache_dir = dir), c("a", "b"))
  expect_equal(length(fs::dir_ls(dir, glob = "*.rds")), 2)
  readr::write_lines(c("z"), f)
  expect_equal(cached_read(f, reader, cache_dir = dir), "z")
  expect_equal(cache_clear(dir), 3)
})

test_that("plots return ggplot objects on fixture data and NULL on empty input", {
  coh <- fixture_cohort()
  fac <- read_cohort_facets(coh, drop_diploid = FALSE)
  qc <- facets_qc_gate(fac$qc)
  snv <- snv_events(read_cohort_maf(coh))
  cnv <- cnv_events(fac$gene, qc)
  sv <- sv_events(read_cohort_sv(coh))
  svg <- sv_gene_events(sv)
  ev <- bind_events(snv, cnv, svg)
  summ <- event_summary(list(recurrence_table(snv, 3, "SNV"), recurrence_table(cnv, 3, "CNV"),
    recurrence_table(svg, 3, "SV")), 3)
  samples <- sample_summary(coh, snv = snv, sv = sv, facets_qc = qc,
    sample_data = read_cohort_qc(coh)$sample_data)
  figs <- list(
    plot_oncoprint(ev, top_n = 10), plot_recurrent_genes(summ, 10), plot_sv_burden(sv),
    plot_sv_type_mix(sv), plot_sv_size(sv), plot_cytoband_gainloss(cytoband_recurrence(cnv, 3)),
    plot_arm_heatmap(cnv_arm_events(fac$arm, qc)), plot_purity_ploidy(qc), plot_tmb(samples)
  )
  for (p in figs) expect_s3_class(p, "ggplot")
  expect_null(plot_sv_burden(sv[0, ]))
  expect_null(plot_oncoprint(ev[0, ]))
  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  save_figs(set_names(figs, str_c("f", seq_along(figs))), pdf_file)
  expect_true(fs::file_exists(pdf_file))
})
