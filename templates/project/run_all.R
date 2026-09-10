# run_all.R
#
# Run the standard WGSCohortAnalysis pipeline for this project:
#
#   Rscript run_all.R            # all stages
#   Rscript run_all.R 03 04      # only stages 03 and 04
#
# Reads 00.PARAMS.yml, loads the toolkit named there (or WCA_HOME), then
# sources scripts/01..06 from the toolkit in order.

params <- yaml::read_yaml("00.PARAMS.yml")
toolkit <- Sys.getenv("WCA_HOME", unset = params$toolkit %||% "")
`%||%` <- function(a, b) if (is.null(a)) b else a
if (!nzchar(toolkit) || !file.exists(file.path(toolkit, "load.R"))) {
  stop("Set `toolkit:` in 00.PARAMS.yml or the WCA_HOME environment variable")
}
source(file.path(toolkit, "load.R"))

stages <- sort(list.files(file.path(toolkit, "scripts"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
wanted <- commandArgs(trailingOnly = TRUE)
if (length(wanted) > 0) {
  stages <- stages[substr(basename(stages), 1, 2) %in% wanted]
}
for (stage in stages) {
  source(stage, echo = FALSE)
}
cat(sprintf("\nrun_all.R finished (%d stages)\n", length(stages)))
