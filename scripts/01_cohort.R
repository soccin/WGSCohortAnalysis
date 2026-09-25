# scripts/01_cohort.R
#
# Stage 1: manifests -> cohort table.
# Reads both Tempo manifests, applies the cohort rules in order (exclude
# regex, project include, dedup by Sig), checks TID uniqueness and set
# symmetry, resolves snv_dir / sv_dir independently, and writes
# cache/<run>/01_cohort/cohort.rds and results/<run>/tables/cohort.xlsx.
#
# Run from a project root: Rscript $WCA_HOME/scripts/01_cohort.R

if (!exists("wca_home", mode = "function")) {
  p0 <- yaml::read_yaml("00.PARAMS.yml")
  source(file.path(Sys.getenv("WCA_HOME", unset = p0$toolkit %||% ""), "load.R"))
}
PARAMS <- wca_read_params("00.PARAMS.yml")
STAGE <- "01_cohort"
stage_banner(STAGE, PARAMS)

workbook <- NULL
if (!is.null(PARAMS$cohort$projects_file)) {
  workbook <- read_projects_file(project_path(PARAMS, PARAMS$cohort$projects_file))
}
selection <- select_cohort_projects(PARAMS$cohort$include_projects,
  PARAMS$cohort$include_ucode, workbook)
include_projects <- selection$projects
projects <- selection$workbook
if (!is.null(projects)) wca_msg("  projects file: {nrow(projects)} projects selected")
if (length(include_projects) == 0) {
  wca_msg("  cohort projects: no project filter, every project in the manifests")
} else {
  wca_msg("  cohort projects ({length(include_projects)}): {str_c(include_projects, collapse = ', ')}")
}

read_side <- function(key) {
  file <- PARAMS$manifests[[key]]
  if (is.null(file)) return(NULL)
  m <- read_manifest(project_path(PARAMS, file)) |>
    apply_cohort_rules(exclude_tid_regex = PARAMS$cohort$exclude_tid_regex,
      include_projects = include_projects)
  wca_msg("  {key} manifest: {nrow(m)} pairs after cohort rules")
  m
}
snv_manifest <- read_side("snv")
sv_manifest <- read_side("sv")

cohort <- build_cohort(snv_manifest, sv_manifest, require_both = isTRUE(PARAMS$cohort$require_both))
wca_msg("  cohort: {nrow(cohort)} samples across {n_distinct(cohort$ProjNo)} projects; ",
  "{sum(cohort$has_snv)} with SNV, {sum(cohort$has_sv)} with SV, {sum(cohort$has_facets)} with FACETS")

stage_save(cohort, PARAMS, STAGE, "cohort")
if (!is.null(projects)) stage_save(projects, PARAMS, STAGE, "projects")

tables <- list(Cohort = cohort)
if (!is.null(projects)) tables$Projects <- projects
write_report_xlsx(tables, fs::path(results_dir(PARAMS, "tables"), "cohort.xlsx"))
wca_msg("  wrote results/{PARAMS$run}/tables/cohort.xlsx")
