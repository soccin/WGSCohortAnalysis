# A string with a comma in it is never one number. readr disagrees by
# default: its grouping mark is "," so "90,0" parses as 900. Every reader in
# the toolkit has to be immune to that, and so does a bare read_tsv() in a
# project script, which is why wca_load() sets the locale for the session.

comma_tsv <- function(dir = tempdir()) {
  path <- fs::path(dir, "comma_pairs.tsv")
  write_lines(c("PR\tSR\tCIPOS\tdepth",
                "90,0\t101,22\t-50,50\t35",
                "86,0\t92,35\t-125,0\t40"), path)
  path
}

test_that("type_convert_silent leaves comma-separated pairs alone", {
  x <- tibble(manta_PR = c("90,0", "86,0"), manta_SR = c("101,22", "92,35"),
              depth = c("35", "40"))
  out <- type_convert_silent(x)
  expect_type(out$manta_PR, "character")
  expect_type(out$manta_SR, "character")
  expect_equal(out$manta_PR, c("90,0", "86,0"))
  # Real numbers still convert; the fix is not "give up on type guessing".
  expect_type(out$depth, "integer")
})

test_that("the session locale protects a bare readr call", {
  d <- read_tsv(comma_tsv(), show_col_types = FALSE, progress = FALSE)
  expect_type(d$PR, "character")
  expect_type(d$SR, "character")
  expect_type(d$CIPOS, "character")
  expect_equal(d$SR, c("101,22", "92,35"))
  expect_true(is.numeric(d$depth))
})

test_that("no reader can turn a comma into a grouping mark", {
  # The structural version of the test above: whatever the column name and
  # whatever else is in the file, text containing a comma comes back as text.
  raw <- read_tsv(comma_tsv(), col_types = cols(.default = "c"), progress = FALSE)
  conv <- type_convert_silent(raw)
  has_comma <- map_lgl(raw, \(x) any(str_detect(x, ",")))
  expect_true(all(map_lgl(conv[has_comma], is.character)))
})

test_that("the locale helpers say what they are", {
  expect_equal(wca_locale()$grouping_mark, "")
  expect_equal(getOption("readr.default_locale")$grouping_mark, "")
  expect_equal(getOption("vroom.default_locale")$grouping_mark, "")
})

test_that("guess_parser under the toolkit locale calls a pair text", {
  expect_equal(guess_parser(c("90,0", "86,0"), locale = wca_locale()), "character")
  expect_equal(guess_parser(c("90", "86"), locale = wca_locale()), "double")
})
