# scripts/04_summaries.R
#
# Stage 4: recurrence tables. Per-type gene recurrence, the wide event
# summary (all events and LoF), cytoband and arm recurrence, SV pairs,
# fusions, SV partners, the per-sample summary and the genes-of-interest
# slices.

if (!exists("wca_home", mode = "function")) {
  p0 <- yaml::read_yaml("00.PARAMS.yml")
  source(file.path(Sys.getenv("WCA_HOME", unset = p0$toolkit %||% ""), "load.R"))
}
PARAMS <- wca_read_params("00.PARAMS.yml")
STAGE <- "04_summaries"
stage_banner(STAGE, PARAMS)

cohort <- stage_load(PARAMS, "01_cohort", "cohort")
snv <- stage_load(PARAMS, "03_events", "snv")
cnv_gene <- stage_load(PARAMS, "03_events", "cnv_gene")
cnv_arm <- stage_load(PARAMS, "03_events", "cnv_arm")
facets_qc <- stage_load(PARAMS, "03_events", "facets_qc")
sv <- stage_load(PARAMS, "03_events", "sv")
svg <- stage_load(PARAMS, "03_events", "svg")
den <- stage_load(PARAMS, "03_events", "denominators")
sample_data <- stage_load(PARAMS, "02_read", "sample_data")
qc_status <- stage_load(PARAMS, "02_read", "qc_status")
genome <- PARAMS$genome
n_samples <- den$n_samples
n_cnv_samples <- den$n_cnv_samples
autosomes_only <- isTRUE(PARAMS$cnv$autosomes_only)

recur <- list(
  snv = recurrence_table(snv, n_samples, "SNV"),
  cnv = recurrence_table(cnv_gene, n_cnv_samples, "CNV"),
  sv = recurrence_table(svg, n_samples, "SV")
)
summary_all <- event_summary(recur, n_samples)

recur_lof <- list(
  snv = recurrence_table(filter(snv, is_lof), n_samples, "SNV-LoF"),
  cnv = recurrence_table(filter(cnv_gene, Dir %in% "loss", filter %in% c("PASS", "RESCUE")), n_cnv_samples, "CNV-loss")
)
summary_lof <- event_summary(recur_lof, n_samples)

recur_gain <- recurrence_table(filter(cnv_gene, Dir %in% "gain", filter %in% c("PASS", "RESCUE")), n_cnv_samples, "CNV-gain")

cnv_bands <- cytoband_recurrence(cnv_gene, n_cnv_samples, genome, autosomes_only = autosomes_only)
arms <- arm_recurrence(cnv_arm, n_cnv_samples, min_frac = PARAMS$cnv$arm_min_frac, genome = genome,
  autosomes_only = autosomes_only)
if (autosomes_only && !all(c(cnv_bands$chrom, arms$chrom) %in% autosomes(genome))) {
  wca_abort("Non-autosomal rows reached the CNV summaries with autosomes_only = TRUE")
}
sv_pairs <- sv_pair_recurrence(sv, n_samples)
fusions <- fusion_table(sv, n_samples)
sv_partners <- sv_partner_table(svg)

samples <- sample_summary(cohort, snv = snv, sv = sv, cnv_gene = cnv_gene,
  facets_qc = if (nrow(facets_qc) > 0) facets_qc else NULL,
  sample_data = sample_data, qc_status = qc_status, autosomes_only = autosomes_only, genome = genome)

goi <- PARAMS$report$genes_of_interest
gene_set <- gene_set_table(summary_all, goi)
gene_set_lof <- gene_set_table(summary_lof, goi)

wca_msg("  genes: {nrow(summary_all)} with any event, {nrow(summary_lof)} with LoF; ",
  "{nrow(cnv_bands)} bands; {nrow(sv_pairs)} SV pairs; {nrow(fusions)} fusion pairs")

for (nm in c("recur", "summary_all", "recur_lof", "summary_lof", "recur_gain", "cnv_bands",
  "arms", "sv_pairs", "fusions", "sv_partners", "samples", "gene_set", "gene_set_lof")) {
  stage_save(get(nm), PARAMS, STAGE, nm)
}
