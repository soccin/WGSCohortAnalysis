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

## gene_model (`read_gene_model()`)

Exon and CDS records for named genes, read from a GENCODE GTF the caller
supplies. One row per feature: gene, transcript, transcript_id, feature
(`exon` or `CDS`), chrom, start, end, strand, exon (number in
transcription order), phase, length. `longest_transcript()` picks one
transcript per gene, `breakpoint_context()` places a position in that
transcript, and `fusion_transcript()` returns the retained exons either
side of a junction with the reading frame.

## Junction anchors (`junction_anchor_support()`)

One row per breakend of a called SV: breakend, chrom, pos, then the
longest exact match between the caller's assembled consensus and the
reference around that breakend (`match_len`, `match_start_in_consensus`,
`match_pos`), `gc_in_match`, `at_fraction` and `informative`. A breakend
whose only match is an A/T tract has `gc_in_match` near zero and is not
located by the data, whatever mapping quality the caller reports.

## Alignment evidence (`bam_reads()`, `bam_soft_clips()`)

`bam_reads()` returns one row per primary alignment in a window: read_id,
qname, flag, chrom, pos, ref_end, strand, mapq, cigar, mate_chrom,
mate_pos, isize, seq. `pos` and `ref_end` bracket the aligned part only,
so soft-clipped bases fall outside them.

`bam_soft_clips()` is long, one row per clip: read_id, side
(`left`/`right`), clip_pos (the reference base the clip abuts), clip_len,
clip_seq, gc_frac, at_frac, plus pos, ref_end, mapq, strand and the mate
columns. A cluster of clips at one base with a common clipped sequence is
an insertion the reference does not have.

`bam_binned_counts()` gives chrom, bin_start, bin_end, n, rel (n over the
median bin) for evenly sampled windows, which is the cheap test of whether
a called deletion removed anything. `bam_base_depth()` gives chrom, pos,
depth over a small window, walking the CIGAR so deletions and skips are
not counted as covered.

`read_pair_orientation()` adds self_reverse, mate_reverse, self_is_left,
orientation (`FR`/`RF`/`FF`/`RR`, read from the leftmost mate) and
sv_signature. FR is the deletion or normal orientation, RF a tandem
duplication, FF and RR the two junctions of an inversion; only FR pairs
between two loci rules an inversion of the segment between them out.
`clip_partner_match()` adds match_forward, match_revcomp and match_frac,
the longest exact match between each clipped sequence and a candidate
partner locus in each orientation, which is what separates a deletion
junction from an inversion junction from sequence that is in neither.

`clip_genotype()` returns one row: n_reads, n_intact, n_clipped, one
count per named clip class, and vaf_clipped. No intact read means
homozygous for the insertion, both means heterozygous, no clips means the
reference allele. `ref_at_profile()` gives chrom, pos, at_frac in sliding
reference windows.

## Gene expression (`read_featurecounts_genes()`)

gene, gene_id, count, length, library_size, CPM, FPKM for named genes from
one featureCounts file. The library size is the sum over every gene in the
file. `featurecounts_beside()` finds the file that sits next to a Forte
metafusion call.

## Recurrence tables

`recurrence_table()` returns Gene, n, PCT, Event with the denominator as an
attribute; `event_summary()` widens a list of them to Gene, Nt, PCT,
`<Event>_n`, `<Event>_PCT`. `cytoband_recurrence()`, `arm_recurrence()`,
`sv_pair_recurrence()` and `fusion_table()` also carry the denominator
attribute, which `write_report_xlsx()` prints in the DataDictionary sheet.
`group_recurrence()` counts samples per key (Gene, PairKey, GenePair)
within named sample groups: `<key>`, n, then `n_<group>`, `PCT_<group>`,
`Samples_<group>` per group; its denominator attribute is the named vector
of group sizes.

Denominators: `n_samples` (cohort size) for SNV and SV; `n_cnv_samples`
(samples passing `facets_qc_gate()`) for CNV.
