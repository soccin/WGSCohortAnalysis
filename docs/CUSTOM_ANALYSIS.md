# Custom analyses with the toolkit

Three worked examples that use the readers and summaries directly rather
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

## Tips

* `read_tempo_maf(file, cols = NULL, row_filter = keep_all_rows())` reads
  everything when the standard column set is not enough.
* `read_tempo_sv(file, detail = "minimal")` is the fast path for counting.
* `sv_support(sv)` adds the cross-caller MAX/MEDIAN read-support columns
  when a confidence filter beyond `NumCallersPass` is wanted.
* `write_sv_bed(sv, "breakpoints.bed")` produces an IGV track of breakends.
* `annotate_cytoband(df, "chrom", "pos")` and `annotate_arm()` work on any
  table with a chromosome and a position.
