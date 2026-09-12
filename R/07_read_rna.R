# 07_read_rna.R
#
# Gene level expression, for the questions a DNA cohort cannot answer on
# its own: is the gene transcribed at all, and is the partner of a
# predicted fusion transcribed. featureCounts output is the input because
# it is what the RNA pipeline already writes beside its fusion calls.

#' Gene level counts for named genes from a featureCounts table
#'
#' The last column of a featureCounts file is the count, whatever the BAM
#' was called, so it is taken by position. The library size is the sum
#' over every gene in the file, not over the genes asked for.
#'
#' @param path featureCounts output (`#` comment line, then a header).
#' @param gene_ids Named character vector or list: names are symbols,
#'   values are the Geneid used in the file.
#' @return Tibble: gene, gene_id, count, length, library_size, CPM, FPKM.
#'   Genes absent from the file come back with a count of zero.
read_featurecounts_genes <- function(path, gene_ids) {
  ids <- unlist(gene_ids)
  if (length(ids) == 0) wca_abort("read_featurecounts_genes(): gene_ids is empty")
  if (!file_exists(path)) wca_abort("featureCounts file not found: {path}")
  d <- read_tsv(path, comment = "#", locale = wca_locale(),
                show_col_types = FALSE, progress = FALSE)
  require_cols(d, c("Geneid", "Length"), "featureCounts table")
  names(d)[ncol(d)] <- "count"
  total <- sum(d$count, na.rm = TRUE)
  tibble(gene = names(ids), gene_id = unname(ids)) |>
    left_join(d |> select(Geneid, count, length = Length), by = c("gene_id" = "Geneid")) |>
    mutate(count = coalesce(count, 0),
           library_size = total,
           CPM = count / total * 1e6,
           FPKM = count / (total / 1e6) / (length / 1e3))
}

#' The featureCounts file that sits beside a Forte per-sample output
#'
#' Forte writes `<sample>/featurecounts/*.gene.featureCounts.txt` next to
#' `<sample>/metafusion/<sample>.final.cff`, so the fusion call locates
#' the counts without a second manifest.
#'
#' @param fusion_file Path to a metafusion `.final.cff`.
#' @return Path, or NA when nothing matches.
featurecounts_beside <- function(fusion_file) {
  dir <- fs::path(fs::path_dir(fs::path_dir(fusion_file)), "featurecounts")
  if (!fs::dir_exists(dir)) return(NA_character_)
  hit <- fs::dir_ls(dir, glob = "*.gene.featureCounts.txt", fail = FALSE)
  if (length(hit) == 0) NA_character_ else as.character(hit[1])
}
