# scripts/06_report.R
#
# Stage 6: xlsx workbooks.
#   <prefix>_CohortEvents_<yymmdd>.xlsx  Samples, FacetsQC, AllEvents, LossOfFunc,
#       SNVGenes, CNVGenes, CNVGains, CNVBands, CNVArms, SVGenes, SVPairs, Fusions,
#       SVPartners, GeneSet, GeneSetLoF, DataDictionary
#   <prefix>_GeneEvents_<yymmdd>.xlsx    SNV, CNV, SV detail for the genes of interest
# Also copies 00.PARAMS.yml and writes TOOLKIT_VERSION into results/<run>/.

if (!exists("wca_home", mode = "function")) {
  p0 <- yaml::read_yaml("00.PARAMS.yml")
  source(file.path(Sys.getenv("WCA_HOME", unset = p0$toolkit %||% ""), "load.R"))
}
PARAMS <- wca_read_params("00.PARAMS.yml")
STAGE <- "06_report"
stage_banner(STAGE, PARAMS)

s <- function(name) stage_load(PARAMS, "04_summaries", name)
e <- function(name) stage_load(PARAMS, "03_events", name)
den <- e("denominators")
cut <- PARAMS$report$pct_cutoff
goi <- PARAMS$report$genes_of_interest
recur <- s("recur")
summary_all <- s("summary_all")

with_den <- function(x, d) { attr(x, "denominator") <- d; x }

cohort_events <- list(
  Samples = s("samples"),
  FacetsQC = e("facets_qc"),
  AllEvents = summary_all |>
    filter(coalesce(SV_PCT > cut$sv, FALSE) | coalesce(SNV_PCT > cut$snv, FALSE) | coalesce(CNV_PCT > cut$cnv, FALSE)) |>
    with_den(den$n_samples),
  LossOfFunc = s("summary_lof") |> filter(PCT > cut$snv) |> with_den(den$n_samples),
  SNVGenes = recur$snv |> filter(PCT > cut$snv) |> with_den(den$n_samples),
  CNVGenes = recur$cnv |> filter(PCT > cut$cnv) |> with_den(den$n_cnv_samples),
  CNVGains = s("recur_gain") |> filter(PCT > cut$cnv) |> with_den(den$n_cnv_samples),
  CNVBands = s("cnv_bands"),
  CNVArms = s("arms"),
  SVGenes = recur$sv |> filter(PCT > cut$sv) |> with_den(den$n_samples),
  SVPairs = s("sv_pairs"),
  Fusions = s("fusions"),
  SVPartners = s("sv_partners"),
  GeneSet = s("gene_set") |> with_den(den$n_samples),
  GeneSetLoF = s("gene_set_lof") |> with_den(den$n_samples)
)

notes <- c(
  cohort = str_glue("{den$n_samples} samples; CNV denominator {den$n_cnv_samples} FACETS QC-passing samples"),
  cutoffs = str_glue("AllEvents keeps a gene when SNV_PCT > {cut$snv}, SV_PCT > {cut$sv} or CNV_PCT > {cut$cnv}"),
  cnv_rule = str_glue("cn_call preset '{PARAMS$cnv$call_preset}', diploid rule '{PARAMS$cnv$diploid_rule}', |dipLogR| <= {PARAMS$cnv$max_abs_diplogr}"),
  sv_rule = str_glue("SVs with NumCallersPass >= {PARAMS$sv$min_callers_pass}; one row per gene and event"),
  genes_of_interest = str_c(goi, collapse = ", "),
  toolkit = str_glue("WGSCohortAnalysis {toolkit_version()}")
)

stamp <- date_stamp()
out_dir <- results_dir(PARAMS)
events_file <- fs::path(out_dir, str_glue("{PARAMS$prefix}_CohortEvents_{stamp}.xlsx"))
write_report_xlsx(cohort_events, events_file, notes = notes)
wca_msg("  wrote {fs::path_rel(events_file, PARAMS$root)}")

gene_events <- list(
  SNV = e("snv") |> filter(Gene %in% goi) |> arrange(Gene, Sample),
  CNV = e("cnv_gene") |> filter(Gene %in% goi) |> arrange(Gene, Sample),
  SV = e("svg") |> filter(Gene %in% goi) |> arrange(Gene, Sample)
)
gene_file <- fs::path(out_dir, str_glue("{PARAMS$prefix}_GeneEvents_{stamp}.xlsx"))
write_report_xlsx(gene_events, gene_file, notes = notes["genes_of_interest"])
wca_msg("  wrote {fs::path_rel(gene_file, PARAMS$root)}")

record_run_provenance(PARAMS)
wca_msg("  wrote TOOLKIT_VERSION and 00.PARAMS.yml copy")
