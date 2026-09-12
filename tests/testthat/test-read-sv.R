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

test_that("bedpe_breakends_in_window finds a hit from either breakend", {
  dir <- fs::path(tempdir(), "bedpe_scan")
  fs::dir_create(dir)
  hdr <- "#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tTYPE"
  row <- function(ca, sa, cb, sb, id, type) {
    str_c(ca, sa, sa + 1L, cb, sb, sb + 1L, id, ".", "+", "-", type, sep = "\t")
  }
  a <- fs::path(dir, "s1.bedpe")
  write_lines(c(hdr,
    row("6", 51847748L, "6", 101194031L, "ev_A_side", "DEL"),   # window at B
    row("6", 101180470L, "6", 101213829L, "ev_both", "DEL"),    # window at A and B
    row("6", 5000L, "9", 9000L, "ev_far", "BND")), a)
  b <- fs::path(dir, "s2.bedpe")
  write_lines(c(hdr, row("1", 100L, "2", 200L, "ev_none", "BND")), b)

  paths <- c(S1 = a, S2 = b)
  h <- bedpe_breakends_in_window(paths, "6", 101173580L, 101214440L)
  expect_equal(sort(h$ID), c("ev_A_side", "ev_both"))
  expect_true(all(h$Sample == "S1"))
  expect_equal(h$START_B[h$ID == "ev_A_side"], 101194031L)
  expect_type(h$START_A, "integer")
  expect_equal(h$TYPE[h$ID == "ev_both"], "DEL")

  # A sample with nothing in the window contributes no rows, and a window
  # with no hit anywhere returns a typed empty tibble rather than an error.
  expect_false("S2" %in% h$Sample)
  empty <- bedpe_breakends_in_window(paths, "22", 1L, 1000L)
  expect_equal(nrow(empty), 0)
  expect_true(all(c("Sample", "START_A") %in% names(empty)))

  # Missing files are skipped; all missing is an error.
  h2 <- bedpe_breakends_in_window(c(paths, S3 = fs::path(dir, "gone.bedpe")),
                                  "6", 101173580L, 101214440L)
  expect_equal(nrow(h2), nrow(h))
  expect_error(bedpe_breakends_in_window(c(X = fs::path(dir, "gone.bedpe")),
                                         "6", 1L, 2L), "no file")
})
