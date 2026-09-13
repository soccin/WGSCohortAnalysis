mini_gtf <- function() wca_file("tests", "fixtures", "miniModel", "miniModel.gtf")

mini_model <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- read_gene_model(mini_gtf(), c("GA", "GB"))
    cache
  }
})

test_that("read_gene_model keeps only the named genes and types the columns", {
  m <- mini_model()
  expect_setequal(unique(m$gene), c("GA", "GB"))
  expect_false("GAX" %in% m$gene)                      # the decoy must not leak in
  expect_setequal(unique(m$feature), c("exon", "CDS"))
  expect_equal(unique(m$chrom), "1")                   # chr prefix stripped
  expect_type(m$start, "integer")
  expect_type(m$exon, "integer")
  expect_equal(sum(m$feature == "exon"), 10)           # 4 + 2 for GA, 4 for GB
  expect_equal(m |> filter(gene == "GB") |> pull(strand) |> unique(), "-")
  expect_equal(m |> filter(transcript == "GA-201", feature == "exon", exon == 4) |> pull(length), 200L)

  expect_equal(nrow(read_gene_model(mini_gtf(), "GA", features = "CDS")), 6)
  expect_error(read_gene_model(mini_gtf(), "NOSUCHGENE"), "No records")
  expect_error(read_gene_model("/no/such/file.gtf", "GA"), "not found")
  expect_error(read_gene_model(mini_gtf(), character(0)), "at least one")
  expect_error(read_gene_model(mini_gtf(), "GA|GB"), "alphanumeric")
})

test_that("longest_transcript prefers the transcript with the most exons", {
  tx <- longest_transcript(mini_model())
  expect_equal(tx[["GA"]], "GA-201")
  expect_equal(tx[["GB"]], "GB-201")
})

test_that("breakpoint_context reports exons in transcription order on both strands", {
  m <- mini_model()

  intronic <- breakpoint_context(m, "GA-201", 2500)
  expect_equal(intronic$region, "intron")
  expect_equal(intronic$exon, 2L)
  expect_equal(intronic$next_exon, 3L)
  expect_equal(intronic$dist_prev, 401L)               # 2500 - 2099
  expect_equal(intronic$dist_next, 500L)               # 3000 - 2500

  exonic <- breakpoint_context(m, "GA-201", 2050)
  expect_equal(exonic$region, "exon")
  expect_equal(exonic$exon, 2L)
  expect_true(is.na(exonic$next_exon))

  # On the minus strand the previous exon is the one at the higher coordinate.
  minus <- breakpoint_context(m, "GB-201", 7500)
  expect_equal(minus$exon, 2L)
  expect_equal(minus$next_exon, 3L)
  expect_equal(minus$dist_prev, 500L)                  # 8000 - 7500
  expect_equal(minus$dist_next, 401L)                  # 7500 - 7099

  expect_equal(breakpoint_context(m, "GA-201", 500)$region, "upstream")
  expect_equal(breakpoint_context(m, "GA-201", 99999)$region, "downstream")
  expect_equal(breakpoint_context(m, "GB-201", 99999)$region, "upstream")
  expect_error(breakpoint_context(m, "NOPE-201", 100), "no exons")
})

test_that("retained_exons keeps the promoter side for 5p on either strand", {
  m <- mini_model()
  expect_equal(retained_exons(m, "GA-201", 2500, "5p")$exon, c(1L, 2L))
  expect_equal(retained_exons(m, "GA-201", 2500, "3p")$exon, c(3L, 4L))
  expect_equal(retained_exons(m, "GB-201", 7500, "5p")$exon, c(1L, 2L))
  expect_equal(retained_exons(m, "GB-201", 7500, "3p")$exon, c(3L, 4L))
})

test_that("fusion_transcript computes retained exons and the reading frame", {
  m <- mini_model()

  # GA CDS through exon 2 is 50 + 100 = 150 bases, a whole number of codons,
  # and GB exon 3 starts a codon (phase 0), so the junction is in frame.
  f <- fusion_transcript(m, "GA-201", "GB-201", 2500, 7500)
  s <- f$summary
  expect_equal(s$gene5, "GA")
  expect_equal(s$gene3, "GB")
  expect_equal(s$last_exon5, 2L)
  expect_equal(s$first_exon3, 3L)
  expect_equal(s$n_exon5, 2L)
  expect_equal(s$n_exon3, 2L)
  expect_equal(s$cds5, 150)
  expect_equal(s$cds3, 250)
  expect_equal(s$aa5, 50)
  expect_equal(s$aa_total, 133)
  expect_equal(s$frame5, 0)
  expect_equal(s$phase3, 0L)
  expect_true(s$in_frame)
  expect_equal(f$exons5$exon, c(1L, 2L))
  expect_equal(f$exons3$exon, c(3L, 4L))

  # Three exons of GA is 250 bases, one past a codon boundary: out of frame.
  out <- fusion_transcript(m, "GA-201", "GB-201", 3500, 7500)$summary
  expect_equal(out$frame5, 1)
  expect_false(out$in_frame)

  expect_error(fusion_transcript(m, "GA-201", "GB-201", 500, 7500), "retains no exons")
})

test_that("fusion_transcript reads a GTF phase as the bases that finish a codon", {
  m <- mini_model()

  # A transcript joined to itself across one intron is its own native splice,
  # so it must test in frame. Across the three introns of each gene the 5'
  # side leaves 2, 0 and 1 bases over and the next exon has phase 1, 0 and 2.
  native <- tribble(
    ~tx,      ~pos, ~frame5, ~phase3,
    "GA-201", 1500, 2,       1L,
    "GA-201", 2500, 0,       0L,
    "GA-201", 3500, 1,       2L,
    "GB-201", 8500, 2,       1L,                        # minus strand
    "GB-201", 7500, 0,       0L,
    "GB-201", 6500, 1,       2L
  )
  pwalk(native, function(tx, pos, frame5, phase3) {
    s <- fusion_transcript(m, tx, tx, pos, pos)$summary
    expect_equal(s$frame5, frame5)
    expect_equal(s$phase3, phase3)
    expect_true(s$in_frame, label = glue("{tx} native splice at {pos}"))
  })

  # One base left over needs a phase 2 exon: GA exons 1-3 :: GB exon 4.
  expect_true(fusion_transcript(m, "GA-201", "GB-201", 3500, 6500)$summary$in_frame)

  # Two bases left over do not match phase 2: GA exon 1 :: GB exon 4.
  # Testing frame5 == phase3 would call this one in frame.
  off <- fusion_transcript(m, "GA-201", "GB-201", 1500, 6500)$summary
  expect_equal(off$frame5, 2)
  expect_equal(off$phase3, 2L)
  expect_false(off$in_frame)
})

test_that("plot_fusion_exons draws a track per parent plus the chimera", {
  m <- mini_model()
  f <- fusion_transcript(m, "GA-201", "GB-201", 2500, 7500)
  p <- plot_fusion_exons(f, m, title = "GA::GB")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nlevels(built$plot$data$track), 3)
  expect_equal(sum(built$plot$data$kept), 2 + 2 + 4)   # kept on each parent, all of the chimera
})
