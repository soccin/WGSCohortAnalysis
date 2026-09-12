# Tempo output on disk

What the readers expect, from the survey of the KEJ runs (2024 to 2026).

## Directory template

```
<run>/                                  e.g. .../Proj_17297_B/IRIS/out/17329
  somatic/<TUMOR>__<NORMAL>/            one directory per pair
    combined_mutations/
      <PAIR>.somatic.final.maf          filtered somatic MAF (read by default)
      <PAIR>.somatic.unfiltered.maf     ~10x larger; not read
    combined_svs/
      <PAIR>.final.bedpe                ensemble somatic SVs (read)
      <PAIR>.final.clustered.bedpe, <PAIR>.unfiltered.bedpe, circos.html,
      _exposures.tsv, _catalogues.pdf   (SV signatures; not read)
    facets/<PAIR>/
      <PAIR>.facets_qc.txt              purity, ploidy, dipLogR, wgd, fga, facets_qc
      <PAIR>_OUT.txt                    run parameters
      facets<version>/                  e.g. facets0.5.14c1000pc5000
        <PAIR>.gene_level.txt
        <PAIR>.arm_level.txt
        <PAIR>_hisens.cncf.txt, _hisens.seg, _purity.cncf.txt, _purity.seg
        <PAIR>.qc.txt                   per-fit QC (two rows: purity and hisens cval)
    meta_data/<PAIR>.sample_data.txt    newer runs: purity, ploidy, WGD_status, MSI,
                                        SBS signatures, HLA, Number_of_Mutations, TMB
    multiqc/<PAIR>.QC_Status.txt        older runs: sample, Status, Reason (tumor + normal)
    brass/ delly/ manta/ svaba/ mutect2/ strelka2/ conpair/ hrdetect/ neoantigen/
    svclone/ lohhla/                    per-caller output; not read
  cohort_level/<cohort>/
    sample_data.txt                     all pairs; alignment_qc.txt, cna_genelevel.txt,
    alignment_qc.txt                    cna_armlevel.txt, mut_somatic.maf, sv_somatic.bedpe
                                        (concatenations; the per-pair files are used instead)
  bams/ runlog/
```

`tempo_pair_files()` resolves every per-pair file from any one path;
`tempo_cohort_files()` finds the cohort-level sample_data and alignment_qc.

## Module presence by run vintage

| Vintage | Example | multiqc/QC_Status | meta_data/sample_data | SV callers |
|---|---|---|---|---|
| 2024 NoSVs/DownSV runs | `WGS/KEJ.01/Proj_14887/NoSVs`, `DownSV` | yes | no (cohort_level only) | DownSV only: delly, manta, svaba, brass |
| 2025 IRIS runs | `Projects/Proj_17297_B/IRIS/out/17329` | no | yes | delly, manta, svaba, brass |
| 2025-26 direct runs | `Projects/Proj_17608/out/17608` | no | yes | same |

`read_cohort_qc()` returns empty tables for whichever file is absent and
`sample_summary()` fills what it can.

## The NoSV / DownSV split

For Proj_14887 and the Umich Proj_16840_C samples Tempo was run twice: a
`NoSVs`/`NoSV` run (MAF + FACETS) and a `DownSV` run (bedpe, downsampled
BAMs). The SNV manifest points at the first, the SV manifest at the second,
and `build_cohort()` keeps `snv_dir` and `sv_dir` separately. FACETS and QC
are always taken from `snv_dir`.

## The 17495_I / 17945_I oddity

Proj_17495_I output lives under `out/17945_I/` (digits transposed in the
Tempo project name), and five pairs are also present under
`results/r_002/tempo/17945_I/`. The manifest lists both paths with the same
`Sig`; `apply_cohort_rules()` keeps the first by `distinct(Sig)`.

## File schemas

### somatic.final.maf

Tab-delimited, no comment line in the runs seen (the reader skips leading
`#` lines anyway). 270 columns in 2024-25 runs, 259 in 2025-26 runs: the
difference is the 11 `neo_*` neoantigen columns. Key columns: Hugo_Symbol,
NCBI_Build (`GRCh37`), Chromosome (1-22, X, unprefixed), Start_Position,
Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2,
Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, HGVSp_Short, t_depth,
t_alt_count, n_depth, SIFT, PolyPhen, oncogenic, Hotspot, t_var_freq, tcn,
lcn, cf, purity, ccf_expected_copies, clonality.

### final.bedpe

`##` header lines (fileformat, FILTER, contig), then a `#CHROM_A` header and
43 columns:

```
#CHROM_A START_A END_A CHROM_B START_B END_B ID QUAL STRAND_A STRAND_B TYPE FILTER
NAME_A REF_A ALT_A NAME_B REF_B ALT_B INFO_A INFO_B FORMAT TUMOR NORMAL TUMOR_ID NORMAL_ID
POTENTIAL_CDNA_CONTAMINATION gene1 transcript1 site1 gene2 transcript2 site2 fusion
Cosmic_Fusion_Counts repName-repClass-repFamily:-site1 repName-repClass-repFamily:-site2
CC_Chr_Band CC_Tumour_Types(Somatic) CC_Cancer_Syndrome CC_Mutation_Type
CC_Translocation_Partner DGv_Name-DGv_VarType-site1 DGv_Name-DGv_VarType-site2
```

`ID` is `TEMPO_<CLASS>_<chrA>_<posA>_<chrB>_<posB>_<strands>` with CLASS in
DEL, DUP, INV, INS, TRA; `TYPE` is DEL, DUP, INV, INS, BND. `INFO_A` holds
`key=value` pairs separated by `;` (Callers, NumCallers, NumCallersPass,
SVTYPE, SVLEN, STRANDS, and `<caller>_<field>` entries). `FORMAT` is a
`:`-separated field list (`GT:delly_GQ:...:manta_PR:manta_SR:svaba_...`)
whose values are in `TUMOR` and `NORMAL`; the set of callers per row varies,
so the pivot yields NA where a caller did not call.

Several of those values are a comma-separated pair rather than one number:
manta `PR` and `SR` are `ref,alt`, delly `CIPOS` and `CIEND` are an interval.
They stay as text, and `sv_vaf()` splits them on the comma. See the rule
about grouping marks in `CLAUDE.md` for why that is not negotiable.

### FACETS

* `gene_level.txt`: sample, gene, chrom (1-22, 23 = X), gene_start,
  gene_end, tsg, seg, median_cnlr_seg, segclust, seg_start, seg_end, cf.em,
  tcn.em, lcn.em, cf, tcn, lcn, seg_length, mcn, genes_on_seg, gene_snps,
  gene_het_snps, spans_segs, cn_state, filter.
* `arm_level.txt`: sample, arm (`1p`, `23q`), tcn, lcn, cn_length,
  arm_length, frac_of_arm, cn_state.
* `facets_qc.txt`: tumor_sample_id, path, fit_name, purity_run_* and
  hisens_run_* parameters, is_best_fit, purity, ploidy, dipLogR,
  dipLogR_flag, n_alternative_dipLogR, wgd, fga, n_dip_bal_segs, ...,
  n_amps, n_homdels, ..., facets_qc.
* `_hisens.cncf.txt`: ID, chrom, loc.start, loc.end, seg, num.mark, nhet,
  cnlr.median, mafR, segclust, cnlr.median.clust, mafR.clust, cf, tcn, lcn,
  cf.em, tcn.em, lcn.em. `_hisens.seg`: ID, chrom, loc.start, loc.end,
  num.mark, seg.mean.

### QC

* `QC_Status.txt`: unnamed first column (sample id), Status, Reason; one
  row for the tumor and one for the normal.
* `sample_data.txt`: sample (pair name), purity, ploidy, WGD_status
  (`True`/`False`), MSI_Total_Sites, MSI_Somatic_Sites, MSIscore,
  SBS*.observed / SBS*.pvalue, HLA-A1..C2, Number_of_Mutations, TMB.
* `alignment_qc.txt` (cohort level): Sample, TotalReads, ..., MedianCoverage.

## Manifests

`TID,NID,ProjNo,PATH,Sig`; PATH is the per-pair MAF (SNV manifest) or bedpe
(SV manifest), Sig the md5 of that file. The KEJ manifests are maintained
in `Analysis/00.Manifest/updateManifest.R`; `scan_tempo_outputs()` builds
equivalent ones from a run directory.
