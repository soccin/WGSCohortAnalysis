# Methods

Every filter, threshold and denominator used by the standard pipeline, with
the rationale and the older-script rules that were retired.

## Cohort

1. Read both manifests (`TID,NID,ProjNo,PATH,Sig`).
2. Drop tumors matching `cohort.exclude_tid_regex` (e.g. `_CL` cell lines,
   `CG-AML`).
3. Keep the selected projects (see Cohort selection below).
4. `distinct(Sig)`: a pair reachable by two paths (the 17495_I `out/` and
   `results/r_002/` copies) is read once. `Sig` is the file md5, so equal
   Sig means identical content.
5. TID must be unique in each manifest. With `require_both: true` the SNV
   and SV tumor sets must be identical; otherwise the union is kept and the
   asymmetry logged (`has_snv`, `has_sv`).
6. `snv_dir` and `sv_dir` are resolved independently because the Umich
   samples have SVs in a `DownSV/` run and MAF + FACETS in `NoSV/` runs.

### Cohort selection

Stage 01 builds one list of `ProjNo` values from two sources and takes
their union:

- `include_projects`, as written;
- when `projects_file` is set, the workbook rows whose `Ucode` is in
  `include_ucode`. An empty `include_ucode` selects every row.

If the union is empty, no project filter is applied and every project in
the manifests is kept. `include_ucode` is read only when `projects_file` is
set.

What each setting selects, checked on the 260924 manifests (21 projects,
212 tumors; APTL is 5 projects, 22 tumors):

| `projects_file` | `include_ucode` | `include_projects` | Selects | Tumors |
|---|---|---|---|---|
| null | `[]` | `[]` | every project in the manifests | 212 |
| null | `[APTL]` | `[]` | every project: the Ucode is ignored | 212 |
| null | `[]` | APTL list | the list | 22 |
| workbook | `[]` | `[]` | every project in the workbook | 212 |
| workbook | `[]` | APTL list | every project in the workbook: the list does not narrow it | 212 |
| workbook | `[APTL]` | `[]` | the projects with that Ucode | 22 |
| workbook | `[APTL]` | APTL list | the Ucode projects plus the list | 22 |
| workbook | `[APTX]` (typo) | `[]` | every project in the manifests | 212 |

None of the wrong rows stops with an error. Choose exactly one way:

- **By ProjNo (the default).** `include_projects: [the list]`,
  `projects_file: null`, `include_ucode: []`. Nothing else is consulted.
- **By Ucode.** `projects_file: <workbook>`, `include_ucode: [the code]`,
  `include_projects: []`. First check that every cohort `ProjNo` in the
  manifests has that `Ucode` in the workbook: a project missing from the
  workbook, or listed under another code, is dropped without an error.
  `Ucode` is set per project; a suffixed project does not inherit the base
  project's code.

Never set `projects_file` with `include_ucode` empty.

After stage 01, check the result. The stage prints
`cohort: N samples across P projects`; P must be the number of projects
in the cohort, and the `ProjNo` column of `results/<run>/tables/cohort.xlsx`
must hold exactly those projects.

## SNV

* Reader: columns by name (`data/maf_standard_cols.txt`), all read as
  character then type-converted; `NCBI_Build` must be one of the genome's
  aliases (`GRCh37`, `hg19`, `b37`).
* Non-silent filter (`non_silent()`): drop `Intron, IGR, 3'Flank, 3'UTR,
  5'Flank, 5'UTR, Silent, RNA`, except 5'Flank rows of the genes in
  `snv.keep_5p_flank` (default TERT: the promoter hotspot has no HGVSp and
  would otherwise be lost). This is the one canonical rule; the six older
  variants (`!is.na(HGVSp)`, `!grepl("=$", HGVSp_Short)`, Wurzburg's
  unfiltered read, etc.) are retired.
* `snv_events()` additionally drops rows without `HGVSp_Short` unless they
  are 5'Flank, and applies `snv.min_t_alt_count`, `min_t_depth`, `min_vaf`
  (all 0 by default, matching prior reports).
* LoF classes: Frame_Shift_Del/Ins, Nonsense_Mutation, Splice_Site,
  Translation_Start_Site, Nonstop_Mutation. (AllTCells used only
  `Frame_Shift|Nonsense`; splice and start-loss are now included.)
* Recurrence counts distinct Gene-Sample pairs; a gene hit twice in one
  sample counts once. Denominator `n_samples`.

## CNV

* FACETS gene-level and arm-level tables come from the versioned
  `facets/<PAIR>/facets*/` directory found recursively, so the discovery
  idioms (`dirname(dirname(maf))/facets`, nested `dir_ls`, zip, tar) are
  gone.
* Sex chromosomes are excluded from every CNV table, summary and figure
  by default (`cnv.autosomes_only: true`). FACETS does not model the
  sample's sex, so a normal male X (one copy) is called a chromosome-wide
  loss and an X analysed against the wrong sex is called a gain. Those
  calls are not somatic and, because X carries roughly 800 genes, they
  swamp the recurrence tables. The rule is applied once by
  `filter_autosomes()` and every CNV function (`cnv_events()`,
  `cnv_arm_events()`, `cytoband_recurrence()`, `arm_recurrence()`,
  `sample_summary()`, `plot_cytoband_gainloss()`, `plot_arm_heatmap()`)
  applies it by default, so calling them directly on raw FACETS reads is
  safe. Stages 03 and 04 abort if a non-autosome reaches a CNV table
  while the flag is on, and stage 03 warns when it is off. Turn it off
  only for a deliberate sex-chromosome analysis.
* Chromosome X is coded 23 in FACETS and Y is not emitted. Readers map 23
  to X on load. The autosome filter is `chrom %in% autosomes()` ("1".."22").
  Older scripts that filtered `!chrom %in% c("X","Y")` on raw FACETS output
  silently let 23 through; this is the chrom-23 fix noted in the
  AllTCells validation.
* Sample gate (`facets_qc_gate()`): exclude `abs(dipLogR) > cnv.max_abs_diplogr`
  (1.5; unreliable ploidy). `cnv.require_facets_qc: true` also excludes
  `facets_qc == FALSE`. The default is false to match prior reports, where
  three Umich samples with `facets_qc FALSE` remained in the CNV
  denominator.
* Denominator `n_cnv_samples` = cohort samples passing the gate.
* Diploid rule (`cnv.diploid_rule`): `state` drops `cn_state == "DIPLOID"`
  and keeps CNLOH (default); `tcn` drops `tcn == 2` (the 260719 and
  AllTCells rule, which also dropped CNLOH).
* Direction (`cn_call(preset =)`):
  * `facets` (default): loss = `HOMDEL|HETLOSS|^LOSS|DOUBLE LOSS`,
    gain = `^AMP|^GAIN|GAIN$`, cnloh = `CNLOH`; loss wins, then gain. From
    the AML OGM comparison (`06_driver_genes.R`).
  * `strict`: gain = `GAIN`, loss = `HETLOSS` or `HOMDEL`, everything else
    NA. The 260719 cytoband rule; AMP was not counted there.
  * `tcn`: tcn < 2 loss, tcn > 2 gain, tcn == 2 and lcn == 0 cnloh
    (AllTCells LossOfFunc / CNVAmp rule).
  * `focal`: the `tcn` rule restricted to PASS/RESCUE calls on segments with
    fewer than 10 genes (`recurrentCNV.R` HOMDEL/AMP-style focal rule).
  * Retired: `cn_state %in% c("HOMDEL","AMP")` only (recurrentCNV.R);
    `tcn != 2` as a direction (CellLines).
* LossOfFunc CNV rows: `Dir == "loss"` and `filter` in PASS/RESCUE.
  CNVGains: `Dir == "gain"` and PASS/RESCUE.
* Cytoband recurrence: each gene mapped to its band by `gene_start` in
  `(bandStart, bandEnd]` (UCSC cytoBand is 0-based half-open); a sample
  counts once per band and direction; ranked by `N_gain + N_loss`.
* Arm recurrence: an arm is altered in a sample when the FACETS arm call
  has `Dir` and `frac_of_arm >= cnv.arm_min_frac` (0.5).

## SV

* Reader: all `##` header lines skipped; `INFO_A` and FORMAT pivoted;
  `UUID = TUMOR_ID_ID`; `sv_class` from the ID prefix so BND reads as TRA.
  VAFs per caller: Delly span `DV/(DV+DR)` and junction `RV/(RV+RR)`, Manta
  split-read `ALT/(REF+ALT)`, SvABA `AD/DP`.
* `sv.min_callers_pass` filters on `NumCallersPass`; 0 (default) keeps every
  call as in prior reports, 2 is a medium-confidence set.
* Dedup rule (from 260719): after splitting `GenePair` on `::`, keep one row
  per Gene and UUID, preferring the highest `NumCallersPass`. AllTCells
  260416 did not dedup intragenic events, which double-counted them; that is
  the SV dedup difference noted in its validation.
* Recurrence counts distinct Gene-Sample pairs. Denominator `n_samples`,
  which includes samples with zero SV calls.
* `PairKey` is the sorted `A::B` key (from Wurzburg `ordered_cc`), so
  `sv_pair_recurrence()` is order independent. `GenePair` keeps Tempo's
  breakend order.
* Fusion classes from the `fusion` text: in-frame, out-of-frame,
  protein-fusion, antisense, transcript, intragenic, none, other.

## Report cutoffs

`report.pct_cutoff` = SNV 0.10, SV 0.10, CNV 0.40. AllEvents keeps a gene
when any one type passes its own bar (not the combined PCT); CNV needs the
higher bar because FACETS calls are broad. CNVBands, CNVArms, SVPairs,
Fusions and SVPartners are unfiltered. `genes_of_interest` are matched
exactly, never as a substring (`DAD1` must not pick up `ADAD1`).

## Sample summary

`purity`, `ploidy`, `fga` are blanked when `facets_qc` is FALSE (the
Proj_16840_C report pattern). TMB, Number_of_Mutations and WGD_status come
from Tempo `sample_data.txt` when the run vintage provides it.
