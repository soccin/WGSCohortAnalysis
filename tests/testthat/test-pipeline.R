test_that("run_all.R on a project scaffolded from the template completes end to end", {
  skip_if(Sys.getenv("WCA_SKIP_E2E") == "1")
  home <- wca_home()
  proj <- fs::path(withr::local_tempdir(), "miniProj")
  out <- system2("Rscript", c(fs::path(home, "bin", "wcaNewProject.R"), proj, "mini"),
    stdout = TRUE, stderr = TRUE)
  expect_true(fs::file_exists(fs::path(proj, "00.PARAMS.yml")), info = str_c(out, collapse = "\n"))

  sc <- scan_tempo_outputs(fixture_root())
  paths <- manifest_from_scan(sc, fs::path(proj, "data", "raw"), stamp = "test")
  params <- wca_read_yaml(fs::path(proj, "00.PARAMS.yml"))
  params$manifests <- list(snv = "data/raw/tempoSNVManifest_test.csv", sv = "data/raw/tempoSVManifest_test.csv")
  params$report$genes_of_interest <- list("TP53", "GATA3", "IKZF2")
  params$report$top_n_oncoprint <- 15
  # write_yaml() writes TRUE as `yes`, which wca_read_params() rejects.
  yaml::write_yaml(params, fs::path(proj, "00.PARAMS.yml"),
    handlers = list(logical = yaml::verbatim_logical))

  log <- withr::with_dir(proj, system2("Rscript", "run_all.R", stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status") %||% 0
  expect_equal(status, 0, info = str_c(tail(log, 40), collapse = "\n"))

  res <- fs::path(proj, "results", "run01")
  expect_true(length(fs::dir_ls(res, glob = "*_CohortEvents_*.xlsx")) == 1)
  expect_true(length(fs::dir_ls(res, glob = "*_GeneEvents_*.xlsx")) == 1)
  expect_true(length(fs::dir_ls(fs::path(res, "figures"), glob = "*.pdf")) == 1)
  expect_true(fs::file_exists(fs::path(res, "TOOLKIT_VERSION")))
  expect_true(fs::file_exists(fs::path(res, "00.PARAMS.yml")))
  expect_true(fs::file_exists(fs::path(res, "tables", "cohort.xlsx")))
  expect_gt(length(fs::dir_ls(fs::path(proj, "cache", "reads"), glob = "*.rds")), 0)

  wb <- fs::dir_ls(res, glob = "*_CohortEvents_*.xlsx")
  sheets <- readxl::excel_sheets(wb)
  expect_true(all(c("Samples", "AllEvents", "LossOfFunc", "SNVGenes", "CNVGenes", "CNVBands",
    "CNVArms", "SVGenes", "SVPairs", "Fusions", "SVPartners", "GeneSet", "DataDictionary") %in% sheets))
  samples <- readxl::read_excel(wb, "Samples")
  expect_equal(nrow(samples), 3)

  # CNV outputs never carry sex chromosomes by default
  raw_gene <- readRDS(fs::path(proj, "cache", "run01", "02_read", "facets_gene.rds"))
  expect_true("X" %in% raw_gene$chrom)
  for (nm in c("CNVBands", "CNVArms")) {
    expect_false(any(readxl::read_excel(wb, nm)$chrom %in% c("X", "Y", "23")), info = nm)
  }
  gene_wb <- fs::dir_ls(res, glob = "*_GeneEvents_*.xlsx")
  expect_false(any(readxl::read_excel(gene_wb, "CNV")$chrom %in% c("X", "Y", "23")))
  x_genes <- unique(raw_gene$gene[raw_gene$chrom == "X"])
  for (nm in c("CNVGenes", "CNVGains")) {
    expect_false(any(readxl::read_excel(wb, nm)$Gene %in% x_genes), info = nm)
  }

  # second run of stage 02 must be a cache hit
  t0 <- Sys.time()
  log2 <- withr::with_dir(proj, system2("Rscript", c("run_all.R", "02"), stdout = TRUE, stderr = TRUE))
  expect_equal(attr(log2, "status") %||% 0, 0, info = str_c(tail(log2, 20), collapse = "\n"))
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 60)
})
