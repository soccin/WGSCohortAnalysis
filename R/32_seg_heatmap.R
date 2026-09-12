# 32_seg_heatmap.R
#
# The view IGV gives a .seg file: one row per sample, one coloured
# rectangle per segment, along a stretch of one chromosome. Useful when the
# question is where a segment boundary falls rather than what a gene's copy
# number is, which is what the gene-level tables already answer.
#
# The fill column is the caller's choice. A log2 coverage ratio centred on
# the sample median is the literal IGV view; total copy number minus the
# sample's ploidy is the comparison that survives whole-genome doubling.
# Both are worth drawing, which is why the fill is an argument.

#' Segment heatmap along one chromosome
#'
#' @param segs Tibble with `Sample`, `seg_start`, `seg_end`, the column named
#'   by `fill`, and optionally a logical `loh` column.
#' @param fill Name of the column to colour by, diverging around zero.
#' @param sample_order Sample ids, top row first.
#' @param x_from,x_to Span in base pairs.
#' @param limit Colour scale is clipped at plus and minus this.
#' @param fill_label Label under the colour bar.
#' @param marks Optional tibble `label`, `pos`: a dashed rule and a name.
#' @param spans Optional tibble `label`, `start`, `end`: a two-headed arrow
#'   above the rows, for an interval a caller reported.
#' @param shade Optional tibble `label`, `start`, `end`: a dotted pair, for
#'   the centromere or any other region worth marking but not colouring.
#' @param show_boundaries Draw a tick at every segment edge.
#' @param sample_labels Row labels; defaults to `sample_order`.
#' @param title,subtitle Panel text.
#' @param wrap Subtitle wrap width.
#' @return A ggplot object.
plot_seg_heatmap <- function(segs, fill, sample_order, x_from, x_to,
                             limit = 1.5, fill_label = NULL, marks = NULL,
                             spans = NULL, shade = NULL, show_boundaries = FALSE,
                             sample_labels = NULL, title = NULL, subtitle = NULL,
                             wrap = 120) {
  require_cols(segs, c("Sample", "seg_start", "seg_end", fill), "segs")
  x_from <- as.numeric(x_from)
  x_to <- as.numeric(x_to)
  n <- length(sample_order)
  mb <- function(x) x / 1e6
  in_view <- function(pos) !is.na(pos) & pos >= x_from & pos <= x_to

  dat <- segs |>
    filter(seg_end >= x_from, seg_start <= x_to,
           .data$Sample %in% sample_order) |>
    mutate(
      # Rows run top to bottom in sample order, so y counts down from the top.
      y = n + 1L - match(.data$Sample, sample_order),
      value = pmax(pmin(.data[[fill]], limit), -limit),
      xmin = mb(pmax(seg_start, x_from)),
      xmax = mb(pmin(seg_end, x_to)),
      loh = if ("loh" %in% names(segs)) .data$loh else FALSE
    )
  if (nrow(dat) == 0) return(empty_plot("segment heatmap"))

  has_span <- !is.null(spans) && nrow(spans) > 0
  y_top <- if (has_span) n + 2.0 else n + 1.0

  p <- ggplot(dat) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - 0.42, ymax = y + 0.42,
                  fill = value)) +
    geom_segment(data = filter(dat, loh),
                 aes(x = xmin, xend = xmax, y = y - 0.47, yend = y - 0.47,
                     colour = "minor copy number 0 (loss of heterozygosity)"),
                 linewidth = 1.0)

  if (show_boundaries) {
    edges <- dat |>
      transmute(y, pos = mb(seg_start)) |>
      filter(pos > mb(x_from), pos < mb(x_to))
    p <- p + geom_segment(data = edges,
                          aes(x = pos, xend = pos, y = y - 0.42, yend = y + 0.42,
                              colour = "segment boundary"), linewidth = 0.3)
  }

  # Regions marked but not coloured, drawn as a dotted pair so the segment
  # rectangles are not hidden behind a shaded band.
  if (!is.null(shade) && nrow(shade) > 0) {
    sh <- shade |> filter(in_view(start) | in_view(end))
    if (nrow(sh) > 0) {
      p <- p +
        annotate("segment", x = mb(c(sh$start, sh$end)), xend = mb(c(sh$start, sh$end)),
                 y = 0.3, yend = n + 0.44, colour = "#a8a6a1",
                 linetype = "dotted", linewidth = 0.4) +
        annotate("text", x = mb((sh$start + sh$end) / 2), y = n + 0.54,
                 label = sh$label, size = 2.7, colour = "#8a8884")
    }
  }

  if (!is.null(marks) && nrow(marks) > 0) {
    mk <- marks |> filter(in_view(pos))
    if (nrow(mk) > 0) {
      p <- p +
        annotate("segment", x = mb(mk$pos), xend = mb(mk$pos), y = 0.3, yend = n + 0.66,
                 linetype = "dashed", colour = "#33312e", linewidth = 0.4) +
        annotate("text", x = mb(mk$pos), y = n + 0.86, label = mk$label,
                 size = 3, fontface = "bold", colour = "#33312e")
    }
  }

  if (has_span) {
    p <- p +
      annotate("segment", x = mb(spans$start), xend = mb(spans$end),
               y = n + 1.35, yend = n + 1.35, colour = "#7a2f2f", linewidth = 0.5,
               arrow = arrow(ends = "both", length = unit(2, "mm"), type = "closed")) +
      annotate("text", x = mb((spans$start + spans$end) / 2), y = n + 1.7,
               label = spans$label, size = 3, colour = "#7a2f2f")
  }

  p +
    scale_fill_gradient2(low = "#2166ac", mid = "grey88", high = "#b2182b",
                         midpoint = 0, limits = c(-limit, limit), name = fill_label) +
    scale_colour_manual(values = c("minor copy number 0 (loss of heterozygosity)" = "#000000",
                                   "segment boundary" = "#6b6a67"), name = NULL) +
    scale_y_continuous(breaks = seq_len(n),
                       labels = rev(sample_labels %||% sample_order),
                       limits = c(0.3, y_top), expand = expansion(0)) +
    scale_x_continuous(labels = \(x) str_c(x, " Mb"), expand = expansion(mult = 0.004)) +
    coord_cartesian(xlim = mb(c(x_from, x_to)), clip = "off") +
    guides(fill = guide_colourbar(order = 1),
           colour = guide_legend(order = 2, override.aes = list(linewidth = 1.2))) +
    labs(x = NULL, y = NULL, title = title,
         subtitle = if (is.null(subtitle)) NULL else str_wrap(subtitle, wrap)) +
    theme_wca() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#ededeb", linewidth = 0.3),
          plot.subtitle = element_text(colour = "#52514e"),
          plot.margin = margin(6, 20, 6, 6),
          legend.key.width = unit(12, "mm"), legend.key.height = unit(3.5, "mm"))
}

#' Distance from a position to the nearest segment edge, per sample
#'
#' A real breakpoint at a locus puts a segment boundary within tens of
#' kilobases of it in every sample that carries the event. A called junction
#' whose nearest boundary is megabases away in every sample has no
#' copy-number support, whatever the caller reports.
#'
#' @param segs Tibble with `Sample`, `seg_start`, `seg_end`.
#' @param pos Position in base pairs.
#' @param label Name for the locus.
#' @param window A boundary this close counts as being at the locus.
#' @return Tibble: Locus, Sample, nearest_boundary, distance_bp, distance_mb,
#'   boundary_at_locus.
seg_boundary_distance <- function(segs, pos, label, window = 1e5) {
  require_cols(segs, c("Sample", "seg_start", "seg_end"), "segs")
  segs |>
    group_by(Sample) |>
    reframe(edge = sort(unique(c(seg_start, seg_end)))) |>
    group_by(Sample) |>
    summarize(nearest_boundary = edge[which.min(abs(edge - pos))], .groups = "drop") |>
    mutate(Locus = label,
           distance_bp = as.integer(abs(nearest_boundary - pos)),
           distance_mb = distance_bp / 1e6,
           boundary_at_locus = distance_bp <= window) |>
    select(Locus, Sample, nearest_boundary, distance_bp, distance_mb, boundary_at_locus)
}
