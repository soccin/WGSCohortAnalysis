# 30_plots.R
#
# ggplot2 + patchwork figures. Every plot function returns a ggplot (or
# patchwork) object, or NULL with a message when the input is empty, so the
# figure stage can skip it. Colors are assigned by entity in fixed order and
# never cycled; gain/loss is a two-hue diverging pair.

wca_colors <- list(
  sv_class = c(TRA = "#2a78d6", DEL = "#eb6834", DUP = "#1baf7a", INV = "#eda100", INS = "#e87ba4"),
  cnv = c(gain = "#e34948", loss = "#2a78d6", cnloh = "#4a3aa7"),
  snv = c(SNV = "#008300", `SNV-LoF` = "#0b0b0b"),
  sv = c(SV = "#4a3aa7", Fusion = "#e87ba4"),
  neutral = "#c3c2b7",
  text = "#52514e"
)

#' Toolkit ggplot theme
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
theme_wca <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#e8e8e6", linewidth = 0.3),
      axis.title = ggplot2::element_text(color = wca_colors$text),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

empty_plot <- function(what) {
  wca_msg("  {what}: nothing to plot")
  NULL
}

#' Oncoprint of the long events table
#'
#' Genes are ranked by the number of samples with any event, samples are
#' memo-sorted on the gene order. CNV is the full tile, SNV the inner tile,
#' SV a small mark and fusion a diamond. Frequency bars on top (events per
#' sample) and right (fraction of samples per gene).
#'
#' @param events Tibble from `bind_events()`.
#' @param genes Genes to show; `NULL` picks the top `top_n` by frequency.
#' @param top_n Number of genes when `genes` is `NULL`.
#' @param samples Sample order; `NULL` memo-sorts. Samples absent from
#'   `events` are still drawn (as empty columns) when given here.
#' @param n_samples Denominator for the right bar; defaults to the number
#'   of samples drawn.
#' @return A patchwork object.
plot_oncoprint <- function(events, genes = NULL, top_n = 40, samples = NULL, n_samples = NULL) {
  if (is.null(events) || nrow(events) == 0) return(empty_plot("oncoprint"))
  ev <- events |> filter(!is.na(Gene))
  if (is.null(genes)) {
    genes <- ev |> distinct(Gene, Sample) |> count(Gene, sort = TRUE) |> slice_head(n = top_n) |> pull(Gene)
  }
  ev <- ev |> filter(Gene %in% genes)
  if (nrow(ev) == 0) return(empty_plot("oncoprint"))
  gene_order <- ev |> distinct(Gene, Sample) |> count(Gene) |>
    arrange(desc(n), Gene) |> pull(Gene)
  all_samples <- union(samples %||% character(), unique(ev$Sample))
  mat <- ev |> distinct(Gene, Sample) |> mutate(v = 1L) |>
    tidyr::pivot_wider(names_from = Sample, values_from = v, values_fill = 0L)
  mat <- mat[match(gene_order, mat$Gene), , drop = FALSE]
  m <- as.matrix(mat[, -1, drop = FALSE])
  missing <- setdiff(all_samples, colnames(m))
  if (length(missing) > 0) m <- cbind(m, matrix(0L, nrow(m), length(missing), dimnames = list(NULL, missing)))
  if (is.null(samples)) {
    ord <- do.call(order, c(lapply(seq_len(nrow(m)), \(i) -m[i, ]), list(colnames(m))))
    samples <- colnames(m)[ord]
  }
  n_samples <- n_samples %||% length(samples)

  layer <- function(types) ev |> filter(EventType %in% types) |>
    mutate(Gene = factor(Gene, levels = rev(gene_order)), Sample = factor(Sample, levels = samples))
  cnv <- layer(c("CNV-gain", "CNV-loss", "CNV-cnloh")) |>
    mutate(Alteration = str_remove(EventType, "^CNV-")) |>
    distinct(Gene, Sample, Alteration)
  snv <- layer(c("SNV", "SNV-LoF")) |> group_by(Gene, Sample) |>
    summarize(Alteration = if (any(EventType == "SNV-LoF")) "SNV-LoF" else "SNV", .groups = "drop")
  sv <- layer("SV") |> distinct(Gene, Sample) |> mutate(Type = "SV")
  fus <- layer("Fusion") |> distinct(Gene, Sample) |> mutate(Type = "Fusion")
  grid <- tidyr::expand_grid(Gene = factor(gene_order, levels = rev(gene_order)),
    Sample = factor(samples, levels = samples))
  fills <- c(wca_colors$cnv, wca_colors$snv)
  shapes <- c(SV = 15, Fusion = 18)
  shape_cols <- c(SV = wca_colors$sv[["SV"]], Fusion = wca_colors$sv[["Fusion"]])
  marks <- bind_rows(sv, fus)

  main <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = grid, ggplot2::aes(Sample, Gene), fill = "#f1f1ef", color = "white", linewidth = 0.4) +
    ggplot2::geom_tile(data = cnv, ggplot2::aes(Sample, Gene, fill = Alteration), color = "white", linewidth = 0.4) +
    ggplot2::geom_tile(data = snv, ggplot2::aes(Sample, Gene, fill = Alteration), width = 0.9, height = 0.45) +
    ggplot2::scale_fill_manual(values = fills, breaks = names(fills), name = NULL, drop = TRUE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    theme_wca() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      panel.grid = ggplot2::element_blank(), legend.position = "bottom") +
    ggplot2::labs(x = NULL, y = NULL)
  if (nrow(marks) > 0) {
    main <- main +
      ggplot2::geom_point(data = marks, ggplot2::aes(Sample, Gene, shape = Type, color = Type), size = 2) +
      ggplot2::scale_shape_manual(values = shapes, breaks = names(shapes), name = NULL) +
      ggplot2::scale_color_manual(values = shape_cols, breaks = names(shape_cols), name = NULL)
  }

  top <- ev |> distinct(Gene, Sample, EventType) |> count(Sample) |>
    mutate(Sample = factor(Sample, levels = samples)) |>
    ggplot2::ggplot(ggplot2::aes(Sample, n)) +
    ggplot2::geom_col(fill = wca_colors$neutral, width = 0.8) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    theme_wca() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = "events")
  right <- ev |> distinct(Gene, Sample) |> count(Gene) |>
    mutate(Gene = factor(Gene, levels = rev(gene_order)), PCT = n / n_samples) |>
    ggplot2::ggplot(ggplot2::aes(PCT, Gene)) +
    ggplot2::geom_col(fill = wca_colors$neutral, width = 0.8) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    theme_wca() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "samples", y = NULL)

  patchwork::wrap_plots(top, patchwork::plot_spacer(), main, right,
    ncol = 2, widths = c(6, 1), heights = c(1, 6)) +
    patchwork::plot_annotation(title = str_glue("Oncoprint: top {length(gene_order)} genes, {length(samples)} samples"))
}

#' SV count per sample by class
#'
#' @param sv Tibble from `sv_events()`.
#' @param sample_order Optional sample order; defaults to descending total.
#' @return A ggplot.
plot_sv_burden <- function(sv, sample_order = NULL) {
  if (is.null(sv) || nrow(sv) == 0) return(empty_plot("sv_burden"))
  d <- sv |> distinct(Sample, UUID, sv_class) |> count(Sample, sv_class)
  order <- sample_order %||% (d |> count(Sample, wt = n) |> arrange(desc(n)) |> pull(Sample))
  d |>
    mutate(Sample = factor(Sample, levels = order), sv_class = factor(sv_class, levels = names(wca_colors$sv_class))) |>
    ggplot2::ggplot(ggplot2::aes(Sample, n, fill = sv_class)) +
    ggplot2::geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = wca_colors$sv_class, name = "SV class", drop = FALSE) +
    theme_wca() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    ggplot2::labs(title = "Structural variants per sample", x = NULL, y = "SV count")
}

#' SV class proportions per sample
#'
#' @inheritParams plot_sv_burden
#' @return A ggplot.
plot_sv_type_mix <- function(sv, sample_order = NULL) {
  if (is.null(sv) || nrow(sv) == 0) return(empty_plot("sv_type_mix"))
  d <- sv |> distinct(Sample, UUID, sv_class) |> count(Sample, sv_class)
  order <- sample_order %||% (d |> count(Sample, wt = n) |> arrange(desc(n)) |> pull(Sample))
  d |>
    mutate(Sample = factor(Sample, levels = order), sv_class = factor(sv_class, levels = names(wca_colors$sv_class))) |>
    ggplot2::ggplot(ggplot2::aes(Sample, n, fill = sv_class)) +
    ggplot2::geom_col(position = "fill", width = 0.8, color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = wca_colors$sv_class, name = "SV class", drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    theme_wca() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    ggplot2::labs(title = "SV class mix per sample", x = NULL, y = "fraction of SVs")
}

#' Size distribution of intrachromosomal SVs
#'
#' @param sv Tibble from `sv_events()`.
#' @return A ggplot (log10 span, one panel per class).
plot_sv_size <- function(sv) {
  if (is.null(sv) || nrow(sv) == 0) return(empty_plot("sv_size"))
  d <- sv |> filter(!is_inter, !is.na(span), span > 0) |> distinct(UUID, sv_class, span)
  if (nrow(d) == 0) return(empty_plot("sv_size"))
  d |>
    mutate(sv_class = factor(sv_class, levels = names(wca_colors$sv_class))) |>
    ggplot2::ggplot(ggplot2::aes(span, fill = sv_class)) +
    ggplot2::geom_histogram(bins = 40, color = "white", linewidth = 0.2) +
    ggplot2::scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ggplot2::scale_fill_manual(values = wca_colors$sv_class, name = "SV class", drop = FALSE) +
    ggplot2::facet_wrap(~sv_class, ncol = 1, scales = "free_y") +
    theme_wca() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "SV size (intrachromosomal)", x = "span (bp)", y = "SVs")
}

#' Cytoband gain/loss frequency in genome order
#'
#' Gains are drawn upward, losses downward, one bar per band with
#' chromosome boundaries marked.
#'
#' @param cnv_bands Tibble from `cytoband_recurrence()`.
#' @param genome Genome name.
#' @param autosomes_only Draw chromosomes 1 to 22 only (default; see
#'   `filter_autosomes()`).
#' @return A ggplot.
plot_cytoband_gainloss <- function(cnv_bands, genome = "hg19", autosomes_only = TRUE) {
  if (is.null(cnv_bands) || nrow(cnv_bands) == 0) return(empty_plot("cytoband_gainloss"))
  if (autosomes_only) cnv_bands <- filter_autosomes(cnv_bands, "chrom", genome)
  if (nrow(cnv_bands) == 0) return(empty_plot("cytoband_gainloss"))
  g <- genome_info(genome)
  offsets <- tibble(chrom = names(g$chrom_len), len = unname(g$chrom_len)) |>
    mutate(offset = cumsum(dplyr::lag(len, default = 0)), mid = offset + len / 2)
  d <- cnv_bands |>
    left_join(g$cytobands |> select(Band, bandStart, bandEnd), by = "Band") |>
    left_join(offsets |> select(chrom, offset), by = "chrom") |>
    mutate(x0 = offset + bandStart, x1 = offset + bandEnd) |>
    filter(!is.na(x0))
  if (nrow(d) == 0) return(empty_plot("cytoband_gainloss"))
  long <- bind_rows(
    d |> transmute(x0, x1, Dir = "gain", y = PCT_gain),
    d |> transmute(x0, x1, Dir = "loss", y = -PCT_loss)
  ) |> filter(y != 0)
  used <- offsets |> filter(chrom %in% d$chrom)
  ggplot2::ggplot(long) +
    ggplot2::geom_vline(xintercept = used$offset, color = "#e8e8e6", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, color = wca_colors$text, linewidth = 0.3) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x0, xmax = x1, ymin = 0, ymax = y, fill = Dir)) +
    ggplot2::scale_fill_manual(values = wca_colors$cnv, name = NULL) +
    ggplot2::scale_x_continuous(breaks = used$mid, labels = used$chrom, expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(labels = \(x) scales::percent(abs(x))) +
    theme_wca() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Cytoband gain (up) and loss (down) frequency", x = NULL, y = "fraction of CNV samples")
}

#' Arm-level heatmap (sample x arm)
#'
#' @param cnv_arm Tibble from `cnv_arm_events()`.
#' @param genome Genome name.
#' @param min_frac Minimum `frac_of_arm` for a call to be colored.
#' @param autosomes_only Draw chromosomes 1 to 22 only (default; see
#'   `filter_autosomes()`).
#' @return A ggplot.
plot_arm_heatmap <- function(cnv_arm, genome = "hg19", min_frac = 0.5, autosomes_only = TRUE) {
  if (is.null(cnv_arm) || nrow(cnv_arm) == 0) return(empty_plot("arm_heatmap"))
  if (autosomes_only) cnv_arm <- filter_autosomes(cnv_arm, "chrom", genome)
  if (nrow(cnv_arm) == 0) return(empty_plot("arm_heatmap"))
  arms <- chrom_arms(genome) |> filter(chrom %in% unique(cnv_arm$chrom)) |> pull(arm)
  d <- cnv_arm |>
    filter(arm %in% arms) |>
    mutate(Dir = if_else(coalesce(frac_of_arm >= min_frac, FALSE), Dir, NA_character_)) |>
    distinct(Sample, arm, Dir) |>
    mutate(arm = factor(arm, levels = arms))
  sample_order <- d |> filter(!is.na(Dir)) |> count(Sample) |> arrange(desc(n)) |> pull(Sample)
  sample_order <- c(sample_order, setdiff(unique(d$Sample), sample_order))
  d |>
    mutate(Sample = factor(Sample, levels = rev(sample_order))) |>
    ggplot2::ggplot(ggplot2::aes(arm, Sample, fill = Dir)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = wca_colors$cnv, na.value = "#f1f1ef", name = NULL, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    theme_wca() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 6), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = str_glue("Arm-level calls (frac_of_arm >= {min_frac})"), x = NULL, y = NULL)
}

#' Most recurrently altered genes by event type
#'
#' @param summary Tibble from `event_summary()`.
#' @param top_n Genes to show.
#' @return A ggplot with one bar group per gene.
plot_recurrent_genes <- function(summary, top_n = 30) {
  if (is.null(summary) || nrow(summary) == 0) return(empty_plot("recurrent_genes"))
  pct_cols <- setdiff(str_subset(names(summary), "_PCT$"), "PCT")
  d <- summary |> slice_head(n = top_n) |>
    select(Gene, all_of(pct_cols)) |>
    tidyr::pivot_longer(-Gene, names_to = "Event", values_to = "PCT") |>
    mutate(Event = str_remove(Event, "_PCT$"), Gene = factor(Gene, levels = rev(summary$Gene[seq_len(min(top_n, nrow(summary)))])))
  pal <- c(SNV = wca_colors$snv[["SNV"]], CNV = wca_colors$cnv[["gain"]], SV = wca_colors$sv[["SV"]],
    `SNV-LoF` = wca_colors$snv[["SNV-LoF"]], `CNV-loss` = wca_colors$cnv[["loss"]], `CNV-gain` = wca_colors$cnv[["gain"]])
  ggplot2::ggplot(d, ggplot2::aes(PCT, Gene, fill = Event)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    theme_wca() +
    ggplot2::labs(title = str_glue("Top {nlevels(d$Gene)} recurrently altered genes"), x = "fraction of samples", y = NULL)
}

#' FACETS purity versus ploidy
#'
#' @param facets_qc Gated QC table.
#' @return A ggplot.
plot_purity_ploidy <- function(facets_qc) {
  if (is.null(facets_qc) || nrow(facets_qc) == 0) return(empty_plot("purity_ploidy"))
  d <- facets_qc |>
    mutate(status = case_when(
      "EXCLUDE" %in% names(facets_qc) & coalesce(EXCLUDE, FALSE) ~ "excluded",
      !coalesce(facets_qc, TRUE) ~ "facets_qc FALSE",
      TRUE ~ "pass"))
  wgd <- if ("wgd" %in% names(d)) d$wgd else FALSE
  ggplot2::ggplot(d, ggplot2::aes(ploidy, purity, color = status, shape = coalesce(as.logical(wgd), FALSE))) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::scale_color_manual(values = c(pass = "#2a78d6", `facets_qc FALSE` = "#eda100", excluded = "#e34948"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "WGD") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    theme_wca() +
    ggplot2::labs(title = "FACETS purity and ploidy", x = "ploidy", y = "purity")
}

#' Tumor mutation burden per sample
#'
#' Uses `TMB` from Tempo sample_data when present, else the non-silent
#' mutation count.
#'
#' @param samples Tibble from `sample_summary()`.
#' @return A ggplot.
plot_tmb <- function(samples) {
  if (is.null(samples) || nrow(samples) == 0) return(empty_plot("tmb"))
  if ("TMB" %in% names(samples) && any(!is.na(samples$TMB))) {
    d <- samples |> transmute(Sample, y = TMB); ylab <- "TMB (mutations / Mb)"
  } else if ("n_snv" %in% names(samples)) {
    d <- samples |> transmute(Sample, y = n_snv); ylab <- "non-silent mutations"
  } else {
    return(empty_plot("tmb"))
  }
  d <- d |> filter(!is.na(y)) |> arrange(desc(y)) |> mutate(Sample = factor(Sample, levels = Sample))
  if (nrow(d) == 0) return(empty_plot("tmb"))
  ggplot2::ggplot(d, ggplot2::aes(Sample, y)) +
    ggplot2::geom_col(fill = "#2a78d6", width = 0.8) +
    theme_wca() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    ggplot2::labs(title = "Mutation burden", x = NULL, y = ylab)
}

#' Save one figure as PDF and PNG
#'
#' @param plot ggplot or patchwork.
#' @param path Output path without extension.
#' @param width,height Inches.
#' @param formats Extensions to write.
#' @return Paths written, invisibly.
save_fig <- function(plot, path, width = 10, height = 7, formats = c("pdf", "png")) {
  out <- purrr::map_chr(formats, \(ext) {
    f <- str_c(path, ".", ext)
    ggplot2::ggsave(f, plot, width = width, height = height, dpi = 150, bg = "white")
    f
  })
  invisible(out)
}

#' Save a list of figures as one multi-page PDF plus PNGs
#'
#' @param figs Named list of plots.
#' @param pdf_file Multi-page PDF path.
#' @param png_dir Directory for `<prefix>_<name>.png`; `NULL` skips PNGs.
#' @param prefix PNG file prefix.
#' @param width,height Page size in inches.
#' @return `pdf_file`, invisibly.
save_figs <- function(figs, pdf_file, png_dir = NULL, prefix = "wca", width = 11, height = 8) {
  figs <- purrr::compact(figs)
  if (length(figs) == 0) return(invisible(NULL))
  grDevices::pdf(pdf_file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in figs) print(p)
  if (!is.null(png_dir)) {
    purrr::iwalk(figs, \(p, name) save_fig(p, fs::path(png_dir, str_c(prefix, "_", name)),
      width = width, height = height, formats = "png"))
  }
  invisible(pdf_file)
}
