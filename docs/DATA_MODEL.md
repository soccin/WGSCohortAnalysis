# Data model

The tidy tables every function reads and writes. Sample identity is always
the tumor id as `Sample`, with `NID` (matched normal) alongside. Chromosomes
are unprefixed character, FACETS 23/24 mapped to X/Y, ordered
`1..22, X, Y, MT` by `chrom_factor()`.

## cohort (`build_cohort()`)

One row per tumor/normal pair.

| Column | Type | Meaning |
|---|---|---|
| Sample | chr | tumor id (manifest TID) |
| NID | chr | normal id |
| ProjNo | chr | project number |
| snv_file, sv_file | chr | manifest PATH for the MAF / bedpe (NA when absent) |
| snv_dir, sv_dir | chr | per-pair Tempo directory (two levels above the file), resolved independently |
| has_snv, has_sv, has_facets | lgl | file presence; `has_facets` checks for a gene_level file under `snv_dir` |
| Sig_snv, Sig_sv | chr | manifest md5 |

## maf (`read_tempo_maf()`, `read_cohort_maf()`)

One row per non-silent mutation (default `row_filter = non_silent()`), the
columns of `data/maf_standard_cols.txt` selected by name, typed, plus
`UUID = Tumor_Sample_Barcode:Chromosome:Start_Position:Ref:Alt`. The cohort
reader prepends `Sample`, `NID`.

## snv (`snv_events()`)

Analysis view of `maf`: rows with a protein change or a 5'Flank class.

| Column | Meaning |
|---|---|
| Gene, Sample, NID | identity |
| Alteration | HGVSp_Short |
| VClass | Variant_Classification |
| VAF, Depth, AltCount | t_var_freq, t_depth, t_alt_count |
| Chrom, Pos, Ref, Alt | locus |
| HGVSc, Consequence, Hotspot, oncogenic, clonality, SIFT, PolyPhen | annotations (NA when the schema lacks them) |
| UUID | from `maf` |
| is_lof | VClass in `maf_lof_classes()` |

## sv (`read_tempo_sv()`, `sv_events()`)

One row per SV event. The reader pivots `INFO_A` (`Callers`, `NumCallers`,
`NumCallersPass`, `SVLEN`, caller-specific keys) and `FORMAT`/`TUMOR`/`NORMAL`
(`t_<caller>_<field>`, `n_<caller>_<field>`) into columns, adds
`UUID = TUMOR_ID_ID`, `sv_class` (TRA/DEL/DUP/INV/INS from the ID; `TYPE`
keeps bedpe BND), `GenePair = gene1::gene2`, `is_inter`, `span`, and the VAF
columns from `sv_vaf()`. `sv_events()` adds `PairKey` (sorted pair),
`fusion_class`, `is_fusion`, `bandA`, `bandB` and applies
`min_callers_pass`. Priority column order: `data/sv_report_cols.txt`.

## svg (`sv_gene_events()`)

Gene by SV event: `GenePair` split on `::` so fusions appear under both
partners; one row per Gene and UUID keeping the highest `NumCallersPass`
(so an intragenic event yields one row). Columns: Gene, Sample, NID,
sv_class, fusion, fusion_class, gene1, gene2, site1, site2, CHROM_A, START_A,
CHROM_B, START_B, bandA, bandB, NumCallers, NumCallersPass, Callers, UUID.

## cnv_gene (`cnv_events()`)

Gene by sample, non-diploid, FACETS QC-passing samples, chromosomes 1 to
22 only by default (X and Y FACETS calls are not sex-aware; see
docs/METHODS.md). The raw `read_cohort_facets()$gene` table still holds X.

| Column | Meaning |
|---|---|
| Gene, Sample, NID, chrom, gene_start, gene_end | identity and locus |
| tcn, lcn, mcn, cf | FACETS copy numbers and cellular fraction |
| cn_state, filter, genes_on_seg, spans_segs | FACETS annotations |
| Dir | gain / loss / cnloh / NA from `cn_call(preset =)` |
| focal | PASS/RESCUE and fewer than 10 genes on the segment |

## cnv_arm (`cnv_arm_events()`)

Sample, NID, arm, chrom, tcn, lcn, frac_of_arm, cn_state, Dir.

## facets_qc (`read_facets_qc()`, `facets_qc_gate()`)

Sample, NID, facets_qc, purity, ploidy, dipLogR, wgd, fga, n_segs, the
remaining FACETS QC columns, then `EXCLUDE` and `exclude_reason` from the
gate.

## qc_status, sample_data (`read_cohort_qc()`)

`qc_status`: Sample, NID, id (as written in the file), is_tumor, Status,
Reason; both tumor and normal rows. `sample_data`: Sample, NID, purity,
ploidy, WGD_status, MSIscore, Number_of_Mutations, TMB.

## events (`bind_events()`)

Long table: Sample, Gene, EventType, Detail. EventType in
`SNV, SNV-LoF, CNV-gain, CNV-loss, CNV-cnloh, SV, Fusion`. LoF SNVs appear
under both `SNV` and `SNV-LoF`; fusions under both `SV` and `Fusion`.

## Recurrence tables

`recurrence_table()` returns Gene, n, PCT, Event with the denominator as an
attribute; `event_summary()` widens a list of them to Gene, Nt, PCT,
`<Event>_n`, `<Event>_PCT`. `cytoband_recurrence()`, `arm_recurrence()`,
`sv_pair_recurrence()` and `fusion_table()` also carry the denominator
attribute, which `write_report_xlsx()` prints in the DataDictionary sheet.

Denominators: `n_samples` (cohort size) for SNV and SV; `n_cnv_samples`
(samples passing `facets_qc_gate()`) for CNV.
