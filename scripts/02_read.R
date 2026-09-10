# scripts/02_read.R
#
# Stage 2: per-pair reads through the content-addressed cache/reads/.
# Reads MAF (non-silent rows only), bedpe, FACETS gene/arm/qc, QC_Status
# and sample_data for every cohort pair and stores the raw bound tables in
# cache/<run>/02_read/. A second run is a cache hit and takes seconds.

if (!exists("wca_home", mode = "function")) {
  p0 <- yaml::read_yaml("00.PARAMS.yml")
  source(file.path(Sys.getenv("WCA_HOME", unset = p0$toolkit %||% ""), "load.R"))
}
PARAMS <- wca_read_params("00.PARAMS.yml")
STAGE <- "02_read"
stage_banner(STAGE, PARAMS)

cohort <- stage_load(PARAMS, "01_cohort", "cohort")
reads <- reads_cache_dir(PARAMS)
genome <- PARAMS$genome

wca_msg("Reading MAFs (non-silent; keep 5'Flank for {str_c(PARAMS$snv$keep_5p_flank, collapse = ',')})")
maf <- read_cohort_maf(cohort, cache_dir = reads,
  row_filter = non_silent(keep_5p_flank = PARAMS$snv$keep_5p_flank), genome = genome)
stage_save(maf, PARAMS, STAGE, "maf")

wca_msg("Reading SVs")
sv_raw <- read_cohort_sv(cohort, cache_dir = reads, genome = genome)
stage_save(sv_raw, PARAMS, STAGE, "sv_raw")

wca_msg("Reading FACETS")
facets <- read_cohort_facets(cohort, cache_dir = reads, drop_diploid = FALSE)
stage_save(facets$gene, PARAMS, STAGE, "facets_gene")
stage_save(facets$arm, PARAMS, STAGE, "facets_arm")
stage_save(facets$qc, PARAMS, STAGE, "facets_qc_raw")

wca_msg("Reading QC")
qc <- read_cohort_qc(cohort, cache_dir = reads)
stage_save(qc$qc_status, PARAMS, STAGE, "qc_status")
stage_save(qc$sample_data, PARAMS, STAGE, "sample_data")

missing <- cohort |>
  transmute(Sample,
    maf = !Sample %in% maf$Sample, sv = !Sample %in% sv_raw$Sample,
    facets = !Sample %in% facets$gene$Sample) |>
  filter(maf | sv | facets)
if (nrow(missing) > 0) {
  wca_msg("  NOTE: {nrow(missing)} samples contributed no rows to at least one table",
    " (zero calls or missing file); see 02_read/missing.rds")
}
stage_save(missing, PARAMS, STAGE, "missing")
