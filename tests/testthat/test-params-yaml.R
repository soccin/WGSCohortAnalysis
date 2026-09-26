# YAML 1.1, which the yaml package follows, reads y, n, yes, no, on and off
# as booleans. In a parameters file those letters are biology: Y is
# tyrosine and a chromosome, N is asparagine and an unknown base. A hotspot
# written `ref_aa: Y` came back as TRUE and broke a report. The toolkit
# reads yaml so that only true and false are logical, and stops when a
# logical parameter is written any other way.

yaml_file <- function(lines) {
  path <- withr::local_tempfile(fileext = ".yml", .local_envir = parent.frame())
  write_lines(lines, path)
  path
}

test_that("letters and yes/no words stay text; true and false stay logical", {
  y <- wca_read_yaml(yaml_file(c(
    "aa: [Y, N, y, n]",
    "chrom: Y",
    "words: [yes, no, on, off, Yes, NO]",
    "flags: [true, false, True, FALSE]",
    "bases: [A, C, G, T]"
  )))
  expect_identical(unlist(y$aa), c("Y", "N", "y", "n"))
  expect_identical(y$chrom, "Y")
  expect_identical(unlist(y$words), c("yes", "no", "on", "off", "Yes", "NO"))
  expect_identical(unlist(y$flags), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(unlist(y$bases), c("A", "C", "G", "T"))
})

test_that("the yaml package alone would turn Y and N into booleans", {
  # The reason wca_read_yaml() exists. If this ever fails, the yaml package
  # has moved to YAML 1.2 and the handlers are no longer needed.
  y <- yaml::yaml.load("aa: [Y, N]")
  expect_identical(unlist(y$aa), c(TRUE, FALSE))
})

test_that("wca_read_params keeps amino acids and chromosomes as text", {
  f <- yaml_file(c(
    "run: run01",
    "hotspots:",
    "  - {gene: STAT3, codon: 640, ref_aa: Y}",
    "  - {gene: IDH2, codon: 172, ref_aa: R}",
    "  - {gene: KRAS, codon: 116, ref_aa: N}",
    "locus: {chrom: Y, start: 100}"
  ))
  p <- wca_read_params(f)
  expect_identical(map_chr(p$hotspots, "ref_aa"), c("Y", "R", "N"))
  expect_identical(p$locus$chrom, "Y")
})

test_that("a logical parameter written yes or on stops the read", {
  f <- yaml_file(c("cohort:", "  require_both: yes"))
  err <- expect_error(wca_read_params(f), class = "wca_user_error")
  expect_match(conditionMessage(err), "cohort.require_both must be true or false")
  f2 <- yaml_file(c("cnv:", "  autosomes_only: on"))
  expect_error(wca_read_params(f2), "cnv.autosomes_only", class = "wca_user_error")
})

test_that("true and false still set logical parameters", {
  p <- wca_read_params(yaml_file(c("cohort:", "  require_both: false", "cnv:", "  autosomes_only: TRUE")))
  expect_false(p$cohort$require_both)
  expect_true(p$cnv$autosomes_only)
  # Defaults fill what the file leaves out, and pass the check.
  expect_true(wca_read_params(yaml_file("run: run01"))$cohort$require_both)
})

test_that("the template parameters read cleanly", {
  p <- wca_read_params(fs::path(wca_home(), "templates", "project", "00.PARAMS.yml"))
  expect_true(is.logical(p$cohort$require_both))
  expect_true(is.logical(p$cnv$autosomes_only))
})
