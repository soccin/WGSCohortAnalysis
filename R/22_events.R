# 22_events.R
#
# Turn the raw reads into the tidy analysis tables (see docs/DATA_MODEL.md):
# snv, cnv_gene, cnv_arm, sv, svg and the long events table.

#' SNV analysis table from a cohort MAF
#'
#' Rows without a protein change (`HGVSp_Short` NA) are dropped unless they
#' are promoter (5'Flank) rows, which the non-silent filter only lets
#' through for the keep list (TERT).
#'
#' @param maf Tibble from `read_cohort_maf()` (must have `Sample`).
#' @param min_t_alt_count,min_t_depth,min_vaf Optional evidence thresholds;
#'   0 disables.
#' @return Tibble: Gene, Sample, NID, Alteration, VClass, VAF, Depth,
#'   AltCount, Chrom, Pos, Ref, Alt, HGVSc, Consequence, Hotspot, oncogenic,
#'   clonality, SIFT, PolyPhen, UUID.
snv_events <- function(maf, min_t_alt_count = 0, min_t_depth = 0, min_vaf = 0) {
  if (nrow(maf) == 0) return(snv_empty())
  require_cols(maf, c("Sample", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short",
    "t_depth", "t_alt_count"), "maf")
  if (!"NID" %in% names(maf)) maf$NID <- NA_character_
  if (!"t_var_freq" %in% names(maf)) maf <- mutate(maf, t_var_freq = t_alt_count / t_depth)
  opt <- function(col) if (col %in% names(maf)) maf[[col]] else rep(NA, nrow(maf))
  maf |>
    mutate(
      Hotspot = opt("Hotspot"), oncogenic = opt("oncogenic"), clonality = opt("clonality"),
      SIFT = opt("SIFT"), PolyPhen = opt("PolyPhen"), HGVSc = opt("HGVSc"),
      Consequence = opt("Consequence"), UUID = opt("UUID")
    ) |>
    transmute(
      Gene = Hugo_Symbol, Sample, NID,
      Alteration = HGVSp_Short,
      VClass = Variant_Classification,
      VAF = t_var_freq, Depth = t_depth, AltCount = t_alt_count,
      Chrom = Chromosome, Pos = Start_Position,
      Ref = Reference_Allele, Alt = Tumor_Seq_Allele2,
      HGVSc, Consequence, Hotspot, oncogenic, clonality, SIFT, PolyPhen, UUID
    ) |>
    filter(!is.na(Alteration) | VClass == "5'Flank") |>
    filter(coalesce(AltCount >= min_t_alt_count, TRUE),
      coalesce(Depth >= min_t_depth, TRUE),
      coalesce(VAF >= min_vaf, TRUE)) |>
    mutate(is_lof = VClass %in% maf_lof_classes()) |>
    arrange(Gene, Sample)
}

snv_empty <- function() {
  tibble(Gene = character(), Sample = character(), NID = character(),
    Alteration = character(), VClass = character(), VAF = numeric(),
    Depth = numeric(), AltCount = numeric(), Chrom = character(), Pos = numeric(),
    Ref = character(), Alt = character(), HGVSc = character(), Consequence = character(),
    Hotspot = character(), oncogenic = character(), clonality = character(),
    SIFT = character(), PolyPhen = character(), UUID = character(), is_lof = logical())
}

#' Gene-level CNV table from FACETS gene calls
#'
#' @param facets_gene Tibble from `read_cohort_facets()$gene`.
#' @param facets_qc Gated QC table (`facets_qc_gate()`); samples with
#'   `EXCLUDE` are dropped. `NULL` keeps every sample.
#' @param genome Genome name.
#' @param autosomes_only Keep chromosomes 1 to 22 only (default and
#'   strongly recommended; see `filter_autosomes()` for why X is unreliable
#'   in FACETS output). Set to `FALSE` only for a deliberate sex-chromosome
#'   analysis.
#' @param diploid_rule `"state"` drops `cn_state == "DIPLOID"` (keeps CNLOH);
#'   `"tcn"` drops `tcn == 2` (the 260719 rule); `"none"` keeps all rows.
#' @param preset Passed to `cn_call()`.
#' @return Tibble: Gene, Sample, NID, chrom, gene_start, gene_end, tcn, lcn,
#'   mcn, cf, cn_state, filter, genes_on_seg, spans_segs, Dir, focal.
cnv_events <- function(facets_gene, facets_qc = NULL, genome = "hg19", autosomes_only = TRUE,
                       diploid_rule = c("state", "tcn", "none"), preset = "facets") {
  diploid_rule <- match.arg(diploid_rule)
  if (nrow(facets_gene) == 0) return(cnv_empty())
  require_cols(facets_gene, c("Sample", "gene", "chrom", "gene_start", "tcn", "lcn", "cn_state", "filter"), "facets_gene")
  out <- facets_gene
  if (!is.null(facets_qc)) out <- filter(out, Sample %in% cnv_samples(facets_qc))
  if (autosomes_only) out <- filter_autosomes(out, "chrom", genome)
  out <- switch(diploid_rule,
    state = filter(out, cn_state != "DIPLOID"),
    tcn = filter(out, tcn != 2),
    none = out
  )
  if (!"NID" %in% names(out)) out$NID <- NA_character_
  opt <- function(col) if (col %in% names(out)) out[[col]] else rep(NA, nrow(out))
  out |>
    mutate(gene_end = opt("gene_end"), mcn = opt("mcn"), cf = opt("cf"),
      genes_on_seg = opt("genes_on_seg"), spans_segs = opt("spans_segs")) |>
    transmute(
      Gene = gene, Sample, NID, chrom, gene_start, gene_end,
      tcn, lcn, mcn, cf, cn_state, filter, genes_on_seg, spans_segs,
      Dir = cn_call(cn_state, tcn, lcn, filter, genes_on_seg, preset = preset),
      focal = filter %in% c("PASS", "RESCUE") & coalesce(genes_on_seg < 10, FALSE)
    ) |>
    arrange(Gene, Sample)
}

cnv_empty <- function() {
  tibble(Gene = character(), Sample = character(), NID = character(), chrom = character(),
    gene_start = numeric(), gene_end = numeric(), tcn = numeric(), lcn = numeric(),
    mcn = numeric(), cf = numeric(), cn_state = character(), filter = character(),
    genes_on_seg = numeric(), spans_segs = logical(), Dir = character(), focal = logical())
}

#' Arm-level CNV table from FACETS arm calls
#'
#' @inheritParams cnv_events
#' @param facets_arm Tibble from `read_cohort_facets()$arm`.
#' @return Tibble: Sample, NID, arm, chrom, tcn, lcn, frac_of_arm, cn_state, Dir.
cnv_arm_events <- function(facets_arm, facets_qc = NULL, genome = "hg19", autosomes_only = TRUE,
                           preset = "facets") {
  if (nrow(facets_arm) == 0) {
    return(tibble(Sample = character(), NID = character(), arm = character(), chrom = character(),
      tcn = numeric(), lcn = numeric(), frac_of_arm = numeric(), cn_state = character(), Dir = character()))
  }
  out <- facets_arm
  if (!is.null(facets_qc)) out <- filter(out, Sample %in% cnv_samples(facets_qc))
  if (autosomes_only) out <- filter_autosomes(out, "chrom", genome)
  if (!"NID" %in% names(out)) out$NID <- NA_character_
  out |>
    transmute(Sample, NID, arm, chrom, tcn, lcn, frac_of_arm, cn_state,
      Dir = cn_call(cn_state, tcn, lcn, preset = preset)) |>
    arrange(Sample, chrom_factor(chrom, genome), arm)
}

#' Fusion class from the Tempo `fusion` annotation
#'
#' @param fusion Character vector.
#' @return One of in-frame, out-of-frame, protein-fusion, antisense,
#'   transcript, intragenic, none, other.
fusion_class <- function(fusion) {
  case_when(
    is.na(fusion) | fusion == "-" | fusion == "" ~ "none",
    str_detect(fusion, "Protein Fusion") & str_detect(fusion, "in frame") ~ "in-frame",
    str_detect(fusion, "Protein Fusion") & str_detect(fusion, "out of frame") ~ "out-of-frame",
    str_detect(fusion, "Protein Fusion") ~ "protein-fusion",
    str_detect(fusion, "Antisense") ~ "antisense",
    str_detect(fusion, "Transcript Fusion") ~ "transcript",
    str_detect(fusion, "^(Deletion|Duplication|Inversion)") ~ "intragenic",
    TRUE ~ "other"
  )
}

#' Order-independent gene pair key
#'
#' @param gene1,gene2 Character vectors.
#' @return `"A::B"` with the genes sorted.
gene_pair_key <- function(gene1, gene2) {
  a <- pmin(gene1, gene2)
  b <- pmax(gene1, gene2)
  str_c(a, "::", b)
}

#' SV analysis table
#'
#' @param sv Tibble from `read_cohort_sv()`.
#' @param genome Genome name for cytoband annotation.
#' @param min_callers_pass Keep SVs with `NumCallersPass >= this`; 0 keeps all.
#' @return `sv` with `Sample`, `PairKey`, `fusion_class`, `is_fusion`,
#'   `bandA`, `bandB` added and filtered.
sv_events <- function(sv, genome = "hg19", min_callers_pass = 0) {
  if (nrow(sv) == 0) return(sv)
  require_cols(sv, c("Sample", "GenePair", "gene1", "gene2", "CHROM_A", "START_A",
    "CHROM_B", "START_B", "NumCallersPass", "ID"), "sv")
  if (!"sv_class" %in% names(sv)) sv <- mutate(sv, sv_class = sv_class(ID))
  sv |>
    filter(coalesce(NumCallersPass >= min_callers_pass, min_callers_pass <= 0)) |>
    mutate(
      PairKey = gene_pair_key(gene1, gene2),
      fusion_class = fusion_class(fusion),
      is_fusion = str_detect(coalesce(fusion, "-"), "Fusion")
    ) |>
    annotate_cytoband("CHROM_A", "START_A", genome, name = "bandA") |>
    annotate_cytoband("CHROM_B", "START_B", genome, name = "bandB") |>
    select(GenePair, Sample, sv_class, any_of("TYPE"), is_inter, span, fusion, fusion_class,
      is_fusion, gene1, gene2, site1, site2, CHROM_A, START_A, END_A, CHROM_B, START_B, END_B,
      bandA, bandB, NumCallers, NumCallersPass, Callers, PairKey, matches("VAF$"), UUID,
      everything())
}

#' Gene by SV-event table
#'
#' Splits `GenePair` so a fusion appears under both partners, then keeps one
#' row per Gene and event (`UUID`), preferring the highest `NumCallersPass`.
#' Intragenic events (gene1 == gene2) therefore yield one row.
#'
#' @param sv Tibble from `sv_events()`.
#' @return Tibble: Gene, Sample, NID, sv_class, fusion, fusion_class, gene1,
#'   gene2, site1, site2, CHROM_A, START_A, CHROM_B, START_B, bandA, bandB,
#'   NumCallers, NumCallersPass, Callers, UUID.
sv_gene_events <- function(sv) {
  if (nrow(sv) == 0) {
    return(tibble(Gene = character(), Sample = character(), NID = character(),
      sv_class = character(), fusion = character(), fusion_class = character(),
      gene1 = character(), gene2 = character(), site1 = character(), site2 = character(),
      CHROM_A = character(), START_A = numeric(), CHROM_B = character(), START_B = numeric(),
      bandA = character(), bandB = character(), NumCallers = numeric(),
      NumCallersPass = numeric(), Callers = character(), UUID = character()))
  }
  if (!"NID" %in% names(sv)) sv$NID <- NA_character_
  sv |>
    tidyr::separate_longer_delim(GenePair, delim = "::") |>
    rename(Gene = GenePair) |>
    select(Gene, Sample, NID, sv_class, fusion, fusion_class, gene1, gene2, site1, site2,
      CHROM_A, START_A, CHROM_B, START_B, bandA, bandB,
      NumCallers, NumCallersPass, Callers, UUID) |>
    arrange(desc(NumCallersPass)) |>
    distinct(Gene, UUID, .keep_all = TRUE)
}

#' Long events table across SNV, CNV and SV
#'
#' @param snv,cnv_gene,svg Tables from `snv_events()`, `cnv_events()`,
#'   `sv_gene_events()`; any may be `NULL`.
#' @return Tibble: Sample, Gene, EventType, Detail. EventType is one of SNV,
#'   SNV-LoF, CNV-gain, CNV-loss, CNV-cnloh, SV, Fusion. LoF SNVs and fusions
#'   appear both under their general type and their specific type.
bind_events <- function(snv = NULL, cnv_gene = NULL, svg = NULL) {
  parts <- list()
  if (!is.null(snv) && nrow(snv) > 0) {
    parts$snv <- snv |> transmute(Sample, Gene, EventType = "SNV", Detail = coalesce(Alteration, VClass))
    parts$lof <- snv |> filter(is_lof) |>
      transmute(Sample, Gene, EventType = "SNV-LoF", Detail = coalesce(Alteration, VClass))
  }
  if (!is.null(cnv_gene) && nrow(cnv_gene) > 0) {
    parts$cnv <- cnv_gene |> filter(!is.na(Dir)) |>
      transmute(Sample, Gene, EventType = str_c("CNV-", Dir), Detail = cn_state)
  }
  if (!is.null(svg) && nrow(svg) > 0) {
    parts$sv <- svg |> transmute(Sample, Gene, EventType = "SV", Detail = str_c(sv_class, ":", gene1, "::", gene2))
    parts$fus <- svg |> filter(fusion_class %in% c("in-frame", "out-of-frame", "protein-fusion", "transcript")) |>
      transmute(Sample, Gene, EventType = "Fusion", Detail = fusion)
  }
  if (length(parts) == 0) {
    return(tibble(Sample = character(), Gene = character(), EventType = character(), Detail = character()))
  }
  bind_rows(parts) |> distinct() |> arrange(Gene, Sample, EventType)
}
