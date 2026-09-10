maf_file <- function(id = "CL01") tempo_pair_files(fixture_pair(id))$maf

test_that("read_tempo_maf selects columns by name and keeps only non-silent rows", {
  maf <- read_tempo_maf(maf_file("CL03"))
  expect_true(all(maf_standard_cols() %in% names(maf)))
  expect_true("UUID" %in% names(maf))
  expect_false(any(maf$Variant_Classification %in% c("Silent", "Intron", "IGR")))
  expect_true(any(maf$Hugo_Symbol == "TERT" & maf$Variant_Classification == "5'Flank"))
  expect_true(is.numeric(maf$t_depth))
  expect_true(is.numeric(maf$t_var_freq))
})

test_that("keep_all_rows keeps silent and intronic rows", {
  maf <- read_tempo_maf(maf_file("CL01"), row_filter = keep_all_rows())
  expect_true(any(maf$Variant_Classification == "Silent"))
  expect_true(any(maf$Variant_Classification == "Intron"))
  expect_true(any(maf$Chromosome == "X"))
})

test_that("non_silent drops 5'Flank for genes not in the keep list", {
  maf <- read_tempo_maf(maf_file("CL03"), row_filter = non_silent(keep_5p_flank = character()))
  expect_false(any(maf$Variant_Classification == "5'Flank"))
})

test_that("a missing requested column warns instead of failing", {
  expect_warning(read_tempo_maf(maf_file("CL01"), cols = c(maf_standard_cols(), "not_a_column")),
    "absent")
})

test_that("NCBI_Build is asserted against the genome", {
  expect_error(read_tempo_maf(maf_file("CL01"), genome = "hg38"), "Unknown genome")
  f <- withr::local_tempfile(fileext = ".maf")
  lines <- readr::read_lines(maf_file("CL01"))
  readr::write_lines(str_replace_all(lines, "GRCh37", "GRCh38"), f)
  expect_error(read_tempo_maf(f), "NCBI_Build")
})

test_that("read_cohort_maf binds samples with cohort ids and caches", {
  coh <- fixture_cohort()
  cache <- withr::local_tempdir()
  maf <- read_cohort_maf(coh, cache_dir = cache)
  expect_equal(names(maf)[1:2], c("Sample", "NID"))
  expect_setequal(unique(maf$Sample), coh$Sample)
  expect_equal(length(fs::dir_ls(cache, glob = "*.rds")), 3)
  maf2 <- read_cohort_maf(coh, cache_dir = cache)
  expect_equal(maf, maf2)
})
