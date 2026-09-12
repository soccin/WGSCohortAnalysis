# Custom analyses with the toolkit

Six worked examples that use the readers and summaries directly rather
than the standard pipeline. Each starts from a loaded toolkit:

```r
source(file.path(Sys.getenv("WCA_HOME"), "load.R"))
```

If a standard run exists, prefer `stage_load()` over re-reading:

```r
PARAMS <- wca_read_params("00.PARAMS.yml")
snv <- stage_load(PARAMS, "03_events", "snv")
```

## 1. TP53 per-sample report (from CellLines_260314/report02_tp53.R)

One row per sample saying whether TP53 is altered by SNV, CNV (with LOH) or
SV, plus the detail sheets.

```r
cohort <- build_cohort(read_manifest("data/raw/tempoSNVManifest.csv"),
                       read_manifest("data/raw/tempoSVManifest.csv"))
cache <- "cache/reads"

# Only TP53 rows survive the read, so the cache entry is small and specific.
tp53_filter <- structure(
  function(maf) filter(maf, Hugo_Symbol == "TP53", !Variant_Classification %in% maf_noncoding_classes()),
  wca_filter = list(name = "tp53_only"))
maf <- read_cohort_maf(cohort, cache_dir = cache, row_filter = tp53_filter)
snv <- snv_events(maf)

fac <- read_cohort_facets(cohort, cache_dir = cache, drop_diploid = FALSE)
qc  <- facets_qc_gate(fac$qc)
cnv <- cnv_events(fac$gene |> filter(gene == "TP53"), qc, diploid_rule = "state")

sv  <- read_cohort_sv(cohort, cache_dir = cache) |> sv_events()
svg <- sv_gene_events(sv) |> filter(Gene == "TP53")

summary_tbl <- cohort |>
  select(Sample) |>
  left_join(snv |> group_by(Sample) |> summarize(TP53_Mut = TRUE,
    TP53_Mut_Detail = str_c(Alteration, collapse = ";")), by = "Sample") |>
  left_join(cnv |> transmute(Sample, TP53_CNV = cn_state, TP53_LOH = coalesce(lcn == 0, FALSE)), by = "Sample") |>
  left_join(svg |> group_by(Sample) |> summarize(TP53_SV = TRUE,
    TP53_SV_Detail = str_c(unique(str_c(gene1, "::", gene2)), collapse = ";")), by = "Sample") |>
  mutate(across(c(TP53_Mut, TP53_LOH, TP53_SV), \(x) coalesce(x, FALSE)),
    across(c(TP53_Mut_Detail, TP53_CNV, TP53_SV_Detail), \(x) coalesce(x, "")),
    TP53_Altered = TP53_Mut | TP53_CNV != "" | TP53_LOH | TP53_SV)

write_report_xlsx(list(Summary = summary_tbl, TP53_SNV = snv, TP53_CNV = cnv, TP53_SV = svg),
  str_glue("tp53Report_{date_stamp()}.xlsx"))
```

## 2. Gene-set events matrix (from CellLines_260222/analyzeCellLines.R)

A sample-by-gene matrix for a gene list, with recurrence per gene and the
per-pair SV recurrence the cell-line report produced.

```r
genes <- c("TP53", "CDKN2A", "ARID1A", "PLCG1", "CARD11", "STAT3", "STAT5B", "JAK1", "JAK3")

events <- stage_load(PARAMS, "03_events", "events") |> filter(Gene %in% genes)
n_samples <- stage_load(PARAMS, "03_events", "denominators")$n_samples

matrix_tbl <- events |>
  mutate(EventType = str_remove(EventType, "^SNV-LoF$|^Fusion$") |> na_if("")) |>
  filter(!is.na(EventType)) |>
  distinct(Sample, Gene, EventType) |>
  group_by(Sample, Gene) |>
  summarize(Events = str_c(sort(EventType), collapse = "+"), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Gene, values_from = Events, values_fill = "")

sv <- stage_load(PARAMS, "03_events", "sv")
sv_pairs <- sv_pair_recurrence(sv, n_samples, min_n = 2)
sv_genes <- recurrence_table(stage_load(PARAMS, "03_events", "svg"), n_samples, "SV") |>
  filter(Gene %in% genes)

plot_oncoprint(events, genes = genes)
write_report_xlsx(list(Matrix = matrix_tbl, SVGenes = sv_genes, SVPairs = sv_pairs),
  str_glue("geneSet_{date_stamp()}.xlsx"))
```

## 3. Cohort subset by regex (from recurrentSV_APTLSamps_251011.R)

Recurrent SV pairs restricted to the APTL samples of a larger run, with the
denominator being the subset size.

```r
cohort <- stage_load(PARAMS, "01_cohort", "cohort")
sv <- stage_load(PARAMS, "03_events", "sv")

aptl <- cohort |> filter(str_detect(Sample, "Umich|APTL")) |> pull(Sample)
sv_aptl <- sv |> filter(Sample %in% aptl)

pairs <- sv_pair_recurrence(sv_aptl, n_samples = length(aptl))
per_sample <- sv_aptl |> distinct(Sample, PairKey) |> count(Sample, name = "n_pairs")

write_report_xlsx(list(Samples = per_sample, RecurrentSV = pairs, SV = sv_aptl),
  str_glue("recurrentSV_APTL_{date_stamp()}.xlsx"))
```

The `Samples` sheet and the `(denominator)` line in DataDictionary make it
explicit that PCT is relative to the subset, which the original one-off did
not record.

## 4. Events shared by cell lines and patients (from NKCohort/report01_sharedSV.R)

Genes hit by an SV in every cell line and in at least one patient, and the
full ranking of genes hit in at least one patient when nothing is shared by
all. `group_recurrence()` counts samples per gene within each named group;
the questions are then filters on its output.

```r
cohort <- stage_load(PARAMS, "01_cohort", "cohort")
svg <- stage_load(PARAMS, "03_events", "svg")
sv <- stage_load(PARAMS, "03_events", "sv")

groups <- list(
  CellLine = cohort |> filter(str_detect(Sample, "_CL_D")) |> pull(Sample),
  Patient = cohort |> filter(!str_detect(Sample, "_CL_D")) |> pull(Sample)
)
n_cl <- length(groups$CellLine)

sv_genes <- group_recurrence(svg, groups, detail_cols = "sv_class")
shared <- sv_genes |> filter(n_CellLine == n_cl, n_Patient >= 1)
ranked <- sv_genes |> filter(n_Patient >= 1)          # first row = most frequent

# Same for protein chimeras, keeping Tempo's 5'::3' orientation
chimeras <- sv |>
  filter(fusion_class %in% c("in-frame", "protein-fusion")) |>
  group_recurrence(groups, key_col = "GenePair", detail_cols = c("fusion_class", "fusion"))

write_report_xlsx(list(SVGenes_Shared = shared, SVGenes_Ranked = ranked, Chimeras = chimeras),
  str_glue("sharedSV_{date_stamp()}.xlsx"),
  dictionary = bind_rows(wca_dictionary(), group_recurrence_dictionary(groups)))
```

The per-group sizes travel as the `denominator` attribute, so the
DataDictionary reports `PCT_CellLine = n_CellLine / 3 samples` and so on.

## 5. A fusion resolved to exons (from NKCohort/report02_ascc3.R)

What a junction does to the two genes it joins: which exons stay with each
promoter, whether the reading frame survives, and how big the predicted
chimera is. Exon coordinates come from a GENCODE GTF the project names in
its own parameters, so nothing here depends on an annotation package.

```r
sv <- stage_load(PARAMS, "03_events", "sv")
junction <- sv |> filter(GenePair == "PKHD1::ASCC3") |> slice(1)

model <- read_gene_model(PARAMS$ascc3$gene_model_gtf, c("ASCC3", "PKHD1"))
tx <- longest_transcript(model)

# Where each breakend lands: exon, intron, and how far from either exon.
breakpoint_context(model, tx[["ASCC3"]], junction$START_B)
#> region "intron", exon 10, dist_prev 20410, dist_next 20452

# Both genes are on the minus strand and the junction is +/-, so the
# retained halves are ASCC3 5' and PKHD1 3'.
fus <- fusion_transcript(model, tx[["ASCC3"]], tx[["PKHD1"]],
                         junction$START_B, junction$START_A)
fus$summary        # last_exon5 10, first_exon3 36, in_frame TRUE, aa_total 2736
fus$exons5         # the retained exon rows, for a detail sheet

ggsave("fusion.pdf", plot_fusion_exons(fus, model, "ASCC3::PKHD1"),
       width = 9, height = 4)
```

`fusion_transcript()` decides the frame from CDS lengths and the GENCODE
phase of the 3' partner's first retained coding exon, so it agrees with
Tempo's `fusion_class` without trusting it.

## 6. A called junction checked against the alignments (from NKCohort/report03_ascc3_client.R)

When a recurrent junction looks wrong, the alignments settle it. The
question is whether the reference is missing sequence that the sample
has, which is what a mobile element insertion looks like to a caller that
only models rearrangements.

```r
bam <- "path/to/tumor.bam"

# Reads at the locus, and the parts of them the aligner had to cut away.
reads <- bam_reads(bam, "6", 101193772, 101194305)
clips <- bam_soft_clips(reads, min_len = 15) |>
  filter(clip_pos >= 101194022, clip_pos <= 101194055)
clips |> count(side, clip_pos)
#> left 101194032 (poly-A), right 101194045 ((CCCTCT)n)

# The two clipped sequences are the two ends of an SVA element, and the
# 14 bases between the clip points are its target site duplication.
genotype <- clip_genotype(bam, "6", 101194032, 101194045,
  classes = list(hexamer = \(s) str_detect(s, "CCCTCT"),
                 polyA = \(s) str_count(s, "[AT]") / nchar(s) > 0.9))
genotype   # n_intact 0 means both copies carry it

# A 49 Mb deletion would halve the depth over the interval. This does not.
cov <- bam_binned_counts(bam, "6", 45e6, 108e6, n_bins = 160, bin_width = 2000)
cov |> summarize(inside = median(rel[bin_start > 51.8e6 & bin_end < 101.2e6]))

# The second breakend sits on a tract with no G or C, which is the only
# reference sequence a poly-A read can match.
ref_at_profile(fasta, "6", 51847348, 51848148) |> filter(at_frac == 1)

# Deletion or inversion? Orientation answers the second question: an
# inversion of the segment between two loci makes FF and RR pairs and
# cannot make an FR pair across its own span.
reads |>
  filter(mate_chrom == "6", abs(mate_pos - 51847748) < 5000) |>
  read_pair_orientation() |>
  count(orientation, sv_signature)

# Sequence answers both at once. A deletion needs the clipped tail to be
# the partner locus forwards, an inversion needs it reverse complemented.
# Matching neither rules out both.
clip_partner_match(clips, read_ref_seq(fasta, "6", 51847148, 51848348)) |>
  summarize(clip = median(clip_len), fwd = max(match_forward), rev = max(match_revcomp))
```

Expression answers the other half: `read_featurecounts_genes()` on the
Forte output beside each fusion call shows the gene still transcribed and
the predicted partner not transcribed at all, which is why no RNA fusion
call exists.

## Tips

* `read_tempo_maf(file, cols = NULL, row_filter = keep_all_rows())` reads
  everything when the standard column set is not enough.
* `read_tempo_sv(file, detail = "minimal")` is the fast path for counting.
* `sv_support(sv)` adds the cross-caller MAX/MEDIAN read-support columns
  when a confidence filter beyond `NumCallersPass` is wanted.
* `write_sv_bed(sv, "breakpoints.bed")` produces an IGV track of breakends.
* `annotate_cytoband(df, "chrom", "pos")` and `annotate_arm()` work on any
  table with a chromosome and a position.
