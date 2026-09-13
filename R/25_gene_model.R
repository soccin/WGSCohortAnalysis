# 25_gene_model.R
#
# Exon-level gene models read from a GENCODE GTF, and the breakpoint
# arithmetic that turns a pair of SV breakends into a predicted fusion
# transcript. Nothing here is genome-specific beyond the GTF itself; the
# caller supplies the path, so a project pins its own annotation.
#
# The three entry points are `read_gene_model()` (GTF to a tidy exon
# table), `breakpoint_context()` (which exon or intron a breakend lands
# in) and `fusion_transcript()` (retained exons and reading frame across
# a junction). `plot_fusion_exons()` draws the result.

#' Read exon and CDS records for a set of genes from a GENCODE GTF
#'
#' Reads only the `exon` and `CDS` feature lines and only the named genes.
#' A GENCODE GTF is hundreds of megabytes, so the gene filter runs as a
#' `zcat -f | grep` pipe and R only ever sees the matching lines; that
#' also means the file may be plain or gzipped. Chromosomes are normalized
#' by `normalize_chrom()`.
#'
#' @param gtf Path to a GENCODE GTF (plain or gzipped).
#' @param genes Character vector of gene symbols to keep.
#' @param features Feature types to keep.
#' @return Tibble: gene, transcript, transcript_id, feature, chrom, start,
#'   end, strand, exon, phase, length.
read_gene_model <- function(gtf, genes, features = c("exon", "CDS")) {
  if (!file_exists(gtf)) {
    wca_abort("Gene model GTF not found: {gtf}")
  }
  if (length(genes) == 0) {
    wca_abort("read_gene_model() needs at least one gene symbol")
  }
  if (any(str_detect(genes, "[^A-Za-z0-9._-]"))) {
    wca_abort("Gene symbols must be alphanumeric: {str_c(genes, collapse = ', ')}")
  }
  wanted <- str_glue('gene_name "({str_c(genes, collapse = "|")})"')
  # grep exits non-zero when nothing matches, which is our own error below.
  cmd <- str_glue("zcat -f {shQuote(gtf)} | grep -E {shQuote(wanted)} || true")
  lines <- system(cmd, intern = TRUE)
  if (length(lines) == 0) {
    wca_abort("No records for {str_c(genes, collapse = ', ')} in {gtf}")
  }
  raw <- read_tsv(I(lines), col_names = FALSE, col_types = cols(.default = "c"),
                  progress = FALSE) |>
    select(chrom = 1, feature = 3, start = 4, end = 5, strand = 7, phase = 8, attr = 9) |>
    filter(feature %in% features, str_detect(attr, wanted))

  if (nrow(raw) == 0) {
    wca_abort("No {str_c(features, collapse = '/')} records for {str_c(genes, collapse = ', ')} in {gtf}")
  }

  raw |>
    mutate(
      gene = str_extract(attr, 'gene_name "([^"]+)"', group = 1),
      transcript = str_extract(attr, 'transcript_name "([^"]+)"', group = 1),
      transcript_id = str_extract(attr, 'transcript_id "([^"]+)"', group = 1),
      exon = as.integer(str_extract(attr, "exon_number ([0-9]+)", group = 1)),
      chrom = normalize_chrom(chrom),
      start = as.integer(start),
      end = as.integer(end),
      phase = suppressWarnings(as.integer(phase)),
      length = end - start + 1L
    ) |>
    select(gene, transcript, transcript_id, feature, chrom, start, end, strand,
           exon, phase, length) |>
    arrange(gene, transcript, exon, feature)
}

#' Pick the transcript with the most exons for each gene
#'
#' A deliberately simple canonical choice: GENCODE basic already drops the
#' fragments, and "most exons" then matches the RefSeq transcript Tempo
#' annotates against in every case we have checked. Ties break on the
#' longer CDS, then on the transcript name, so the result is stable.
#'
#' @param model Table from `read_gene_model()`.
#' @return Named character vector, transcript name per gene.
longest_transcript <- function(model) {
  require_cols(model, c("gene", "transcript", "feature", "exon", "length"), "longest_transcript")
  picked <- model |>
    group_by(gene, transcript) |>
    summarize(n_exon = sum(feature == "exon"),
              cds_len = sum(length[feature == "CDS"]),
              .groups = "drop_last") |>
    arrange(desc(n_exon), desc(cds_len), transcript, .by_group = TRUE) |>
    slice(1) |>
    ungroup()
  set_names(picked$transcript, picked$gene)
}

#' Locate a breakpoint within a transcript
#'
#' @param model Table from `read_gene_model()`.
#' @param transcript Transcript name.
#' @param pos Breakpoint position, one based.
#' @return One-row tibble: transcript, gene, strand, region (exon, intron,
#'   upstream or downstream), exon (the exon hit, or the exon before the
#'   breakpoint in transcription order), next_exon, dist_prev, dist_next,
#'   where `prev` and `next` are in transcription order, not coordinate
#'   order.
breakpoint_context <- function(model, transcript, pos) {
  ex <- model |> filter(transcript == .env$transcript, feature == "exon")
  if (nrow(ex) == 0) {
    wca_abort("Transcript {transcript} has no exons in the gene model")
  }
  strand <- ex$strand[1]
  gene <- ex$gene[1]
  inside <- ex |> filter(start <= pos, end >= pos)

  # In transcription order, "before" is lower coordinate on +, higher on -.
  before <- if (strand == "+") ex |> filter(end < pos) else ex |> filter(start > pos)
  after <- if (strand == "+") ex |> filter(start > pos) else ex |> filter(end < pos)
  prev_exon <- if (nrow(before) == 0) NA_integer_ else max(before$exon)
  next_exon <- if (nrow(after) == 0) NA_integer_ else min(after$exon)

  gap <- function(e) {
    if (is.na(e)) return(NA_integer_)
    row <- ex |> filter(exon == e) |> slice(1)
    as.integer(if (pos > row$end) pos - row$end else row$start - pos)
  }

  region <- case_when(
    nrow(inside) > 0 ~ "exon",
    is.na(prev_exon) ~ "upstream",
    is.na(next_exon) ~ "downstream",
    TRUE ~ "intron"
  )

  tibble(
    transcript = transcript, gene = gene, strand = strand, pos = as.integer(pos),
    region = region,
    exon = if (nrow(inside) > 0) inside$exon[1] else prev_exon,
    next_exon = if (nrow(inside) > 0) NA_integer_ else next_exon,
    dist_prev = gap(prev_exon),
    dist_next = gap(next_exon)
  )
}

#' Exons retained on one side of a breakpoint
#'
#' @param model Table from `read_gene_model()`.
#' @param transcript Transcript name.
#' @param pos Breakpoint position.
#' @param side `"5p"` keeps the promoter-proximal exons, `"3p"` keeps the
#'   exons downstream of the breakpoint in transcription order.
#' @return Tibble of the retained `exon` feature rows.
retained_exons <- function(model, transcript, pos, side = c("5p", "3p")) {
  side <- match.arg(side)
  ex <- model |> filter(transcript == .env$transcript, feature == "exon")
  if (nrow(ex) == 0) {
    wca_abort("Transcript {transcript} has no exons in the gene model")
  }
  keep_low <- (ex$strand[1] == "+") == (side == "5p")
  if (keep_low) ex |> filter(end <= pos) else ex |> filter(start >= pos)
}

#' Predicted fusion transcript across an SV junction
#'
#' Takes the two breakends of a junction plus the retained side of each and
#' returns the exon structure and reading frame of the chimera. Frame is
#' computed from CDS lengths: the 5' partner contributes
#' `sum(CDS) %% 3` bases into the junction codon and the 3' partner's first
#' retained coding exon carries a GTF phase. A GTF phase is the number of
#' bases before the first whole codon, the bases that finish a codon left
#' open upstream, so the fusion is in frame when `(3 - frame5) %% 3 ==
#' phase3`. Every native splice of a transcript passes this test.
#'
#' @param model Table from `read_gene_model()`.
#' @param tx5,tx3 Transcript names of the 5' and 3' partners.
#' @param pos5,pos3 Breakpoint position in each partner.
#' @return List with `exons5`, `exons3` (tibbles), and a one-row `summary`
#'   tibble: gene5, gene3, last_exon5, first_exon3, n_exon5, n_exon3,
#'   cds5, cds3, aa5, aa3, aa_total, frame5, phase3, in_frame.
fusion_transcript <- function(model, tx5, tx3, pos5, pos3) {
  ex5 <- retained_exons(model, tx5, pos5, "5p")
  ex3 <- retained_exons(model, tx3, pos3, "3p")
  if (nrow(ex5) == 0 || nrow(ex3) == 0) {
    wca_abort("One side of the junction retains no exons ({tx5} / {tx3})")
  }
  cds <- model |> filter(feature == "CDS")
  cds5 <- cds |> filter(transcript == tx5, exon %in% ex5$exon)
  cds3 <- cds |> filter(transcript == tx3, exon %in% ex3$exon)

  len5 <- sum(cds5$length)
  len3 <- sum(cds3$length)
  frame5 <- len5 %% 3
  first3 <- cds3 |> arrange(exon) |> slice(1)
  phase3 <- if (nrow(first3) == 0) NA_integer_ else first3$phase

  list(
    exons5 = ex5 |> arrange(exon),
    exons3 = ex3 |> arrange(exon),
    summary = tibble(
      gene5 = ex5$gene[1], gene3 = ex3$gene[1],
      tx5 = tx5, tx3 = tx3,
      last_exon5 = max(ex5$exon), first_exon3 = min(ex3$exon),
      n_exon5 = nrow(ex5), n_exon3 = nrow(ex3),
      cds5 = len5, cds3 = len3,
      aa5 = len5 %/% 3, aa3 = len3 %/% 3,
      aa_total = (len5 + len3) %/% 3,
      frame5 = frame5, phase3 = phase3,
      in_frame = !is.na(phase3) && (3 - frame5) %% 3 == phase3
    )
  )
}

#' Exon-level schematic of a fusion
#'
#' Two tracks: each parent transcript with its exons in transcription
#' order, the breakpoint marked, and the retained exons filled. A third
#' track shows the chimera. Exon width is uniform, not to scale, because
#' the point is the exon boundary the junction falls between.
#'
#' @param fusion Result of `fusion_transcript()`.
#' @param model Table from `read_gene_model()`.
#' @param title Plot title.
#' @return A ggplot object.
plot_fusion_exons <- function(fusion, model, title = NULL) {
  s <- fusion$summary
  track <- function(tx, kept, label, side) {
    model |>
      filter(transcript == tx, feature == "exon") |>
      arrange(exon) |>
      transmute(track = as.character(label), side = side, exon, idx = row_number(),
                kept = exon %in% kept, gene = gene)
  }
  parents <- bind_rows(
    track(s$tx5, fusion$exons5$exon, str_glue("{s$gene5} ({s$tx5})"), "5p"),
    track(s$tx3, fusion$exons3$exon, str_glue("{s$gene3} ({s$tx3})"), "3p")
  )
  chimera <- bind_rows(
    fusion$exons5 |> transmute(exon, gene),
    fusion$exons3 |> transmute(exon, gene)
  ) |>
    mutate(track = as.character(str_glue("chimera {s$gene5}::{s$gene3}")),
           side = "chimera", idx = row_number(), kept = TRUE)

  dat <- bind_rows(parents, chimera) |>
    mutate(track = factor(track, levels = unique(c(parents$track, chimera$track))))

  # The junction sits after the last retained exon of the 5' partner and
  # before the first retained exon of the 3' partner.
  brk <- dat |>
    filter(side %in% c("5p", "3p"), kept) |>
    group_by(track, side) |>
    summarize(x = if (first(side) == "5p") max(idx) + 0.5 else min(idx) - 0.5,
              .groups = "drop")

  # A 67-exon track cannot carry 67 legible labels; keep the ends only.
  labelled <- dat |>
    filter(kept) |>
    group_by(track) |>
    mutate(small = n() <= 20) |>
    group_by(track, gene) |>
    filter(small | idx %in% range(idx)) |>
    ungroup()

  ggplot(dat, aes(x = idx, y = 1)) +
    geom_tile(aes(fill = gene, alpha = kept), width = 0.85, height = 0.6, colour = "grey30",
              linewidth = 0.2) +
    geom_text(data = labelled, aes(label = exon), size = 2.1, colour = "white") +
    geom_vline(data = brk, aes(xintercept = x), linetype = "dashed", colour = "red") +
    facet_wrap(~track, ncol = 1, scales = "free_x", strip.position = "left") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
    labs(title = title, x = "exon, transcription order", y = NULL, fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          legend.position = "bottom")
}
