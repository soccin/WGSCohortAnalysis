# scripts/05_figures.R
#
# Stage 5: figures. One multi-page PDF plus a PNG per figure under
# results/<run>/figures/.

if (!exists("wca_home", mode = "function")) {
  p0 <- yaml::read_yaml("00.PARAMS.yml")
  source(file.path(Sys.getenv("WCA_HOME", unset = p0$toolkit %||% ""), "load.R"))
}
PARAMS <- wca_read_params("00.PARAMS.yml")
STAGE <- "05_figures"
stage_banner(STAGE, PARAMS)

events <- stage_load(PARAMS, "03_events", "events")
sv <- stage_load(PARAMS, "03_events", "sv")
cnv_arm <- stage_load(PARAMS, "03_events", "cnv_arm")
facets_qc <- stage_load(PARAMS, "03_events", "facets_qc")
summary_all <- stage_load(PARAMS, "04_summaries", "summary_all")
cnv_bands <- stage_load(PARAMS, "04_summaries", "cnv_bands")
samples <- stage_load(PARAMS, "04_summaries", "samples")
genome <- PARAMS$genome
autosomes_only <- isTRUE(PARAMS$cnv$autosomes_only)

figs <- list(
  oncoprint = plot_oncoprint(events, top_n = PARAMS$report$top_n_oncoprint),
  recurrent_genes = plot_recurrent_genes(summary_all, top_n = 30),
  sv_burden = plot_sv_burden(sv),
  sv_type_mix = plot_sv_type_mix(sv),
  sv_size = plot_sv_size(sv),
  cytoband_gainloss = plot_cytoband_gainloss(cnv_bands, genome = genome, autosomes_only = autosomes_only),
  arm_heatmap = plot_arm_heatmap(cnv_arm, genome = genome, autosomes_only = autosomes_only),
  purity_ploidy = plot_purity_ploidy(facets_qc),
  tmb = plot_tmb(samples)
)
figs <- purrr::compact(figs)

fig_dir <- results_dir(PARAMS, "figures")
pdf_file <- fs::path(fig_dir, str_glue("{PARAMS$prefix}_figures_{date_stamp()}.pdf"))
save_figs(figs, pdf_file = pdf_file, png_dir = fig_dir, prefix = PARAMS$prefix)
wca_msg("  wrote {length(figs)} figures to results/{PARAMS$run}/figures/")
