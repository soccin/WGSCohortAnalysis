# scripts/03_events.R
#
# Stage 3: tidy analysis tables. snv, cnv_gene, cnv_arm, sv, svg, events,
# the gated FACETS QC table and the two denominators n_samples (cohort) and
# n_cnv_samples (FACETS QC pass). Stops if a non-autosome reaches cnv_gene
# when autosomes_only is set.

if (!exists("wca_home", mode = "function")) {
  p0 <- yaml::read_yaml("00.PARAMS.yml")
  source(file.path(Sys.getenv("WCA_HOME", unset = p0$toolkit %||% ""), "load.R"))
}
PARAMS <- wca_read_params("00.PARAMS.yml")
STAGE <- "03_events"
stage_banner(STAGE, PARAMS)

cohort <- stage_load(PARAMS, "01_cohort", "cohort")
maf <- stage_load(PARAMS, "02_read", "maf")
sv_raw <- stage_load(PARAMS, "02_read", "sv_raw")
facets_gene <- stage_load(PARAMS, "02_read", "facets_gene")
facets_arm <- stage_load(PARAMS, "02_read", "facets_arm")
facets_qc_raw <- stage_load(PARAMS, "02_read", "facets_qc_raw")
genome <- PARAMS$genome
cnv_p <- PARAMS$cnv

snv <- snv_events(maf, min_t_alt_count = PARAMS$snv$min_t_alt_count,
  min_t_depth = PARAMS$snv$min_t_depth, min_vaf = PARAMS$snv$min_vaf)
wca_msg("  snv: {nrow(snv)} mutations in {n_distinct(snv$Sample)} samples")

facets_qc <- if (nrow(facets_qc_raw) > 0) {
  facets_qc_gate(facets_qc_raw, max_abs_diplogr = cnv_p$max_abs_diplogr,
    require_facets_qc = isTRUE(cnv_p$require_facets_qc))
} else {
  facets_qc_raw
}
excluded <- if (nrow(facets_qc) > 0) facets_qc |> filter(EXCLUDE) |> pull(Sample) else character()
wca_msg("  FACETS gate: {length(excluded)} excluded ({str_c(excluded, collapse = ', ')})")

diploid_rule <- if (isTRUE(cnv_p$drop_diploid)) cnv_p$diploid_rule else "none"
if (!isTRUE(cnv_p$autosomes_only)) {
  wca_warn("cnv.autosomes_only is FALSE: X/Y FACETS calls are not sex-aware and will inflate CNV counts")
}
cnv_gene <- cnv_events(facets_gene, facets_qc = if (nrow(facets_qc) > 0) facets_qc else NULL,
  genome = genome, autosomes_only = isTRUE(cnv_p$autosomes_only),
  diploid_rule = diploid_rule, preset = cnv_p$call_preset)
if (isTRUE(cnv_p$autosomes_only) && !all(cnv_gene$chrom %in% autosomes(genome))) {
  wca_abort("Non-autosomal rows reached cnv_gene with autosomes_only = TRUE")
}
wca_msg("  cnv_gene: {nrow(cnv_gene)} gene calls in {n_distinct(cnv_gene$Sample)} samples")

cnv_arm <- cnv_arm_events(facets_arm, facets_qc = if (nrow(facets_qc) > 0) facets_qc else NULL,
  genome = genome, autosomes_only = isTRUE(cnv_p$autosomes_only), preset = cnv_p$call_preset)
if (isTRUE(cnv_p$autosomes_only) && !all(cnv_arm$chrom %in% autosomes(genome))) {
  wca_abort("Non-autosomal rows reached cnv_arm with autosomes_only = TRUE")
}

sv <- sv_events(sv_raw, genome = genome, min_callers_pass = PARAMS$sv$min_callers_pass)
svg <- sv_gene_events(sv)
wca_msg("  sv: {nrow(sv)} events in {n_distinct(sv$Sample)} samples; svg: {nrow(svg)} gene rows")

events <- bind_events(snv, cnv_gene, svg)

n_samples <- nrow(cohort)
n_cnv_samples <- if (nrow(facets_qc) > 0) length(intersect(cnv_samples(facets_qc), cohort$Sample)) else 0L
wca_msg("  denominators: n_samples = {n_samples}, n_cnv_samples = {n_cnv_samples}")

stage_save(snv, PARAMS, STAGE, "snv")
stage_save(cnv_gene, PARAMS, STAGE, "cnv_gene")
stage_save(cnv_arm, PARAMS, STAGE, "cnv_arm")
stage_save(facets_qc, PARAMS, STAGE, "facets_qc")
stage_save(sv, PARAMS, STAGE, "sv")
stage_save(svg, PARAMS, STAGE, "svg")
stage_save(events, PARAMS, STAGE, "events")
stage_save(list(n_samples = n_samples, n_cnv_samples = n_cnv_samples), PARAMS, STAGE, "denominators")
