# helper-fixtures.R: shared fixture paths and cohort for the tests.

fixture_root <- function() wca_file("tests", "fixtures", "miniCohort")

fixture_pair <- function(id = "CL01") {
  fs::path(fixture_root(), "out", "miniCohort", "somatic", str_glue("{id}__{id}N"))
}

fixture_cohort <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      sc <- scan_tempo_outputs(fixture_root())
      cache <<- build_cohort(sc$snv, sc$sv)
    }
    cache
  }
})
