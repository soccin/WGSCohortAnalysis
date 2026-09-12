bedpe_file <- function(id = "CL01") tempo_pair_files(fixture_pair(id))$bedpe

test_that("read_tempo_sv parses INFO and FORMAT into columns with numeric VAFs", {
  sv <- read_tempo_sv(bedpe_file("CL01"))
  expect_equal(names(sv)[1:3], c("GenePair", "TUMOR_ID", "TYPE"))
  expect_true(all(c("Callers", "NumCallers", "NumCallersPass", "t_delly_DV", "n_delly_DV",
    "t_manta_SR", "t_svaba_AD") %in% names(sv)))
  expect_true(is.numeric(sv$NumCallersPass))
  expect_true(is.numeric(sv$t_delly_SpanVAF))
  expect_true(is.numeric(sv$t_manta_JuncVAF))
  expect_true(is.numeric(sv$t_svaba_VAF))
  expect_true(all(sv$t_delly_SpanVAF >= 0 & sv$t_delly_SpanVAF <= 1, na.rm = TRUE))
  expect_equal(sv$UUID, str_c(sv$TUMOR_ID, "_", sv$ID))
  expect_true(all(sv$is_inter == (sv$CHROM_A != sv$CHROM_B)))
})

test_that("sv_class comes from the ID and maps BND to TRA", {
  sv <- read_tempo_sv(bedpe_file("CL01"))
  expect_setequal(unique(sv$sv_class), c("TRA", "DEL", "DUP", "INV"))
  expect_true(all(sv$sv_class[sv$TYPE == "BND"] == "TRA"))
  expect_equal(sv_class("TEMPO_INS_1_10_1_20_+-"), "INS")
})

test_that("detail levels nest", {
  full <- read_tempo_sv(bedpe_file("CL02"), detail = "full")
  std <- read_tempo_sv(bedpe_file("CL02"), detail = "standard")
  min <- read_tempo_sv(bedpe_file("CL02"), detail = "minimal")
  expect_true(all(names(min) %in% names(std)))
  expect_true(all(names(std) %in% names(full)))
  expect_lt(ncol(min), ncol(std))
  expect_lt(ncol(std), ncol(full))
})

test_that("an empty bedpe returns NULL", {
  f <- withr::local_tempfile(fileext = ".bedpe")
  lines <- readr::read_lines(bedpe_file("CL01"))
  readr::write_lines(lines[str_starts(lines, "#")], f)
  expect_null(read_tempo_sv(f))
})

test_that("sv_support takes max and median across callers", {
  sv <- read_tempo_sv(bedpe_file("CL01")) |> sv_support()
  expect_true(all(c("max_paired_reads", "median_split_reads", "multi_evidence_type") %in% names(sv)))
  expect_true(all(sv$max_split_reads >= sv$median_split_reads, na.rm = TRUE))
})

test_that("read_cohort_sv binds the three fixtures", {
  sv <- read_cohort_sv(fixture_cohort())
  expect_equal(n_distinct(sv$Sample), 3)
  expect_equal(names(sv)[1:3], c("Sample", "NID", "GenePair"))
})

test_that("comma-separated FORMAT pairs are not parsed as grouped numbers", {
  # manta writes PR and SR as "ref,alt". readr's default grouping mark turns
  # "90,0" into 900 and "101,22" into 10122, silently and only when every
  # other value in the column also parses, so the bug moves with the data.
  x <- tibble(manta_PR = c("90,0", "86,0"), manta_SR = c("101,22", "92,35"),
              depth = c("35", "40"))
  out <- type_convert_silent(x)
  expect_type(out$manta_PR, "character")
  expect_type(out$manta_SR, "character")
  expect_equal(out$manta_PR, c("90,0", "86,0"))
  expect_type(out$depth, "integer")
})
