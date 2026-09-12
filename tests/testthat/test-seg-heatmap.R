mini_segs <- function() {
  tibble(
    Sample = rep(c("S1", "S2"), each = 3),
    seg_start = c(1L, 2000001L, 5000001L, 1L, 3000001L, 6000001L),
    seg_end = c(2000000L, 5000000L, 9000000L, 3000000L, 6000000L, 9000000L),
    lr = c(0, -1.2, 0.4, 0.1, -3.0, 2.5),
    loh = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE)
  )
}

test_that("plot_seg_heatmap draws one row per sample, top row first", {
  p <- plot_seg_heatmap(mini_segs(), "lr", c("S1", "S2"), 1, 9e6)
  expect_s3_class(p, "ggplot")

  b <- ggplot2::ggplot_build(p)
  rects <- b$data[[1]]
  expect_equal(nrow(rects), 6)
  # S1 is first in sample_order, so it is the top row: the larger y.
  ys <- p$data |> distinct(Sample, y)
  expect_equal(ys$y[ys$Sample == "S1"], 2)
  expect_equal(ys$y[ys$Sample == "S2"], 1)
  expect_equal(levels(factor(b$layout$panel_params[[1]]$y$get_labels())),
               levels(factor(c("S1", "S2"))))
})

test_that("the fill is clipped at the limit rather than dropped", {
  # -3.0 and 2.5 are outside a limit of 1.5. A scale limit alone would turn
  # them into NA and leave two segments uncoloured, which reads as missing
  # data rather than as a large change.
  p <- plot_seg_heatmap(mini_segs(), "lr", c("S1", "S2"), 1, 9e6, limit = 1.5)
  expect_equal(min(p$data$value), -1.5)
  expect_equal(max(p$data$value), 1.5)
  expect_false(any(is.na(ggplot2::ggplot_build(p)$data[[1]]$fill)))
})

test_that("segments are clipped to the requested span", {
  p <- plot_seg_heatmap(mini_segs(), "lr", c("S1", "S2"), 2.5e6, 4e6)
  expect_true(all(p$data$xmin >= 2.5))
  expect_true(all(p$data$xmax <= 4.0))
  # Only the segments overlapping the window survive.
  expect_equal(nrow(p$data), 3)
})

test_that("loh, marks, spans, shade and boundaries each add a layer", {
  base <- plot_seg_heatmap(mini_segs() |> mutate(loh = FALSE), "lr",
                           c("S1", "S2"), 1, 9e6)
  n_base <- length(base$layers)

  withmarks <- plot_seg_heatmap(mini_segs(), "lr", c("S1", "S2"), 1, 9e6,
    marks = tibble(label = "GENE", pos = 4e6),
    spans = tibble(label = "called deletion", start = 1e6, end = 8e6),
    shade = tibble(label = "centromere", start = 4.5e6, end = 5.5e6),
    show_boundaries = TRUE)
  expect_gt(length(withmarks$layers), n_base)

  # A mark outside the span is dropped rather than drawn off the panel.
  outside <- plot_seg_heatmap(mini_segs(), "lr", c("S1", "S2"), 1, 9e6,
                              marks = tibble(label = "GENE", pos = 50e6))
  expect_equal(length(outside$layers), n_base)
})

test_that("plot_seg_heatmap needs its columns and survives an empty window", {
  expect_error(plot_seg_heatmap(mini_segs() |> select(-seg_end), "lr",
                                c("S1", "S2"), 1, 9e6), "missing columns")
  expect_error(plot_seg_heatmap(mini_segs(), "nope", c("S1", "S2"), 1, 9e6),
               "missing columns")
  expect_null(plot_seg_heatmap(mini_segs(), "lr", c("S1", "S2"), 20e6, 30e6))
})

test_that("seg_boundary_distance finds the nearest edge per sample", {
  d <- seg_boundary_distance(mini_segs(), 2100000, "locus", window = 2e5)
  expect_equal(nrow(d), 2)
  # S1's segments abut at 2,000,000 / 2,000,001, so the nearest edge is
  # 2,000,001: just under 100 kb away and inside the window.
  s1 <- d |> filter(Sample == "S1")
  expect_equal(s1$nearest_boundary, 2000001)
  expect_equal(s1$distance_bp, 99999L)
  expect_true(s1$boundary_at_locus)
  # S2's nearest edge is 3,000,000, 900 kb away, so it is not at the locus.
  s2 <- d |> filter(Sample == "S2")
  expect_equal(s2$nearest_boundary, 3000000)
  expect_equal(s2$distance_bp, 900000L)
  expect_false(s2$boundary_at_locus)
  expect_equal(s2$distance_mb, 0.9)

  expect_error(seg_boundary_distance(mini_segs() |> select(-seg_start), 1, "x"),
               "missing columns")
})
