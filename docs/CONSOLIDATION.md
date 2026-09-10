# Consolidation map

Every surveyed script and its fate. Paths are relative to the survey copy of
the four trees that were consolidated: `Analysis/`, `Projects/`, `WGS/` and a
laptop tree. Owner directory names are left out on purpose; the originals
stay where they are and nothing was copied into `attic/`.

Fates: **function** (logic lives in the named toolkit function),
**pipeline** (covered by the standard stages), **one-off** (project-specific,
not consolidated; run the original if needed), **out of scope** (listed in
the README), **history** (R history dump, ignore).

## The seed: Analysis/260719 and AllTCells_260416

| Script | Fate |
|---|---|
| `Analysis/260719/report01.R` | pipeline: cohort rules -> `apply_cohort_rules()`, `build_cohort()`; MAF filter -> `non_silent()`; snv -> `snv_events()`; FACETS gate -> `facets_qc_gate()`; autosome rule -> `cnv_events()`; cytobands -> `cytoband_recurrence()`; evts -> `recurrence_table()`; sv_partners -> `sv_partner_table()`; output -> `scripts/06_report.R` |
| `Analysis/260719/tools.R` | `event_summary()`, `write_report_xlsx()` |
| `Analysis/AllTCells_260416/R/readTempoPipelineOutput.R` | `read_tempo_sv()` (bedpe body), `read_cohort_files()` + `cached_read()` (memoised generic reader), `tempo_pair_files()` (ls_sample_output_file) |
| `Analysis/AllTCells_260416/R/reportTools.R` | `sv_vaf()` |
| `Analysis/AllTCells_260416/R/reportCols01` | `data/sv_report_cols.txt` |
| `Analysis/AllTCells_260416/R/utils.R` | `type_convert_silent()` |
| `Analysis/AllTCells_260416/report01.R`, `tools.R` | pipeline (LossOfFunc, CNVAmp -> CNVGains, GeneSet regex -> exact `genes_of_interest`) |
| `Analysis/AllTCells_260416/report01_output.md` | `wca_dictionary()` descriptions |

## Cell lines and AML

| Script | Fate |
|---|---|
| `Analysis/CellLines_260222/analyzeCellLines.R` | `docs/CUSTOM_ANALYSIS.md` example 2; `sv_pair_recurrence()` |
| `Analysis/CellLines_260222/analysisTools.R`, `read_tempo_sv.R` | superseded by `read_cohort_files()`, `read_tempo_sv()` |
| `Analysis/CellLines_260222/attic/report01.R` | `read_sample_data()`; svReportMaster rds -> stage rds |
| `Analysis/CellLines_260314/report01.R` | pipeline |
| `Analysis/CellLines_260314/report02_tp53.R` | `docs/CUSTOM_ANALYSIS.md` example 1 |
| `Analysis/CellLines_260314/scripts/00_build_wgs_inventory.R` | `scan_tempo_outputs()`, `manifest_from_scan()` |
| `Analysis/CellLines_260314/attic/scanWGSResults.R`, `getRNADat.R` | one-off (RNA pairing out of scope) |
| `Analysis/AML_260319/report01.R` | pipeline |
| `Analysis/AML_260319/annotCyto.R` | `annotate_cytoband()` (GRanges version replaced by the join_by idiom) |
| `Analysis/00.Manifest/updateManifest.R` | one-off: KEJ manifest maintenance; `scan_tempo_outputs()` is the toolkit equivalent |

## WGS/ (older KEJ analyses)

| Script | Fate |
|---|---|
| `WGS/recurrentSNVs.R`, `recurrentSNVs_V2.R` | pipeline; `get_sample_data()` -> `read_sample_data()`, `tempo_cohort_files()` |
| `WGS/recurrentCNV.R` | `cn_call(preset = "focal")`, `read_facets_gene()`; find_facets_* -> `tempo_pair_files()` |
| `WGS/recurrentSV.R`, `recurrentSV_APTLSamps_251011.R` | `sv_pair_recurrence()`; `docs/CUSTOM_ANALYSIS.md` example 3 |
| `WGS/getSNVs.R`, `analysis_251002.R` | pipeline |
| `WGS/getDAD1Events.R` | out of scope (DAD1) |
| `WGS/gene_cytoband_hg19.R` | `genome_info()`, `annotate_cytoband()` |
| `WGS/Oncoprints/oncoPrint.R`, `mafOncoPrint.R`, `cnvOncoPrint.R`, `svOncoprint.R` | `plot_oncoprint()` (ggplot2; ComplexHeatmap dropped) |
| `WGS/Oncoprints/lollipop.R` | out of scope |
| `WGS/Oncoprints/cnvGetSEGPass.R`, `read_tempo_sv.R` | `read_facets_seg()`, `read_tempo_sv()` |
| `WGS/SVReport/H1.R`, `SVReport/src/*`, `SVReport/misc/resolve` | pipeline (SV sheets); `read_sample_data()` |
| `WGS/SVReport/DAD1/*` | out of scope (DAD1/BLAT) |
| `WGS/SVReport/Analysis/CellLines_251208` | one-off |
| `WGS/SVComparisons/2025-03-27/*` (venn, recurrentGenes, ogm reader) | out of scope (OGM comparison) |
| `WGS/NormalSVs/*` (compareSVs, reportSV01, getTempoInputs) | one-off (normal-vs-normal SV audit) |
| `WGS/NormalStats/getNormalStats.R` | one-off |
| `WGS/KMSS_Set/genReport_250501.R` | pipeline |
| `WGS/Analysis/CellLinesDel_251222/qcReport.R`, `snvReport.R` | pipeline; coverage QC out of scope |
| `WGS/Analysis/2025-05-15_ABLx/report01.R` | one-off |
| `WGS/UMich/getManifest.R`, `plotMetrics.R` | one-off (project QC) |
| `WGS/ReMap_260130/Set*/fixManifest.R` | one-off |
| `WGS/Inventory/*` | empty in the survey copy |
| `WGS/KEJ.01`, `KEJ.03` Module12 | run directories, not scripts |

## Projects/ (per-project one-offs)

| Script | Fate |
|---|---|
| `Projects/Wurzburg/analysis/01_SNV/tools.R` (`read_maf`) | `read_tempo_maf()` |
| `Projects/Wurzburg/analysis/02_SV/tools.R` (`ordered_cc`, `read_tempo_sv`) | `gene_pair_key()`, `read_tempo_sv()` |
| `Projects/Wurzburg/analysis/03_CopyNumber/*` | out of scope (CNV vs RNA) |
| `Projects/Wurzburg/pass1/*`, `Argos/report01.R` | one-off |
| `Projects/Proj_16840_C/GroupA/NoSV/reportSpecialWGS.R` | `sample_summary()` (blank purity/ploidy on failed QC); portal AF join out of scope |
| `Projects/Proj_16840_C/meta/fixSampleIDs.R`, `Proj_16840_D/meta`, `Proj_17495_E/Map/fixPatientIDs.R` | one-off (id fixes) |
| `Projects/Proj_16840_G/mutationCount.R` | `sample_summary()` |
| `Projects/Proj_17929_G/recurAnalysis.R` | pipeline |
| `Projects/Proj_15673/makeSarekInput.R`, `Proj_15966/*`, `Proj_16840_P/TempoBAM/*`, `Proj_17608/*`, `Proj_17495_I/parse_size.R`, `make_pairing.R`, `Proj_17297_B/Eos-2.3.2/SVReport` | one-off (pipeline setup, uploads, pairing) |

## Laptop tree (Mac)

| Script | Fate |
|---|---|
| `Desktop/Work/KEJ/TempoSVOutput/compute_depth_aggregates.R` | `sv_support()` |
| `Desktop/Work/KEJ/TempoSVOutput/reformatDepthReport.R` | one-off |
| `Dropbox/Work/KEJ/DAD1_Event/2025-08-25/dumpSVBED.R` | `sv_to_bed()`, `write_sv_bed()` |
| `Dropbox/Work/KEJ/DAD1_Event/2025-07-31/addDAD1Dist.R` | out of scope (DAD1) |
| `Dropbox/Work/KEJ/WGS_AML/260623/analysis/R/00_common.R` | `genome_info()` (hard-coded lengths replaced), `chrom_factor()` |
| `.../05_cnv_arm_chromosome.R` | `chrom_arms()`, `arm_recurrence()` (frac_of_arm >= 0.5) |
| `.../06_driver_genes.R` | `cn_call(preset = "facets")` regexes |
| `.../01_cohort_overview.R`, `02_sv_concordance.R`, `03_sv_fusions_clinical.R`, `04_cnv_concordance.R`, `07_figures.R` | out of scope (OGM comparison) |
| `Dropbox/Work/KEJ/WGSStats/getStats.R` | one-off (coverage stats) |
| `Dropbox/Work/KEJ/WGS_Facets/geneTable.R` | `read_facets_gene()` |
| `Dropbox/Work/KEJ/WGS_Facets/fout__default/.../_Rhistory_*.R` | history |
| `Dropbox/Work/KEJ/2024-03-28/plotDellyCalls.R` | out of scope (Delly comparison) |
| `Dropbox/Work/KEJ/2025-10-07/analysis_251007.R` | one-off |
| `Desktop/Work/KEJ/WGS/260625/loadSegs.R`, `gene_FBXO45_CNV.R` | `read_facets_seg()`, custom example pattern |
| `Desktop/Work/KEJ/WGS/260625/make_invoice.R` | out of scope (invoice) |
| `Desktop/Work/KEJ/AML_WGS/DataRepo/extractImpactHeme.R`, `reformatClinData.R` | out of scope (IMPACT, clinical) |

## Retired rules (see METHODS.md)

| Old rule | Replacement |
|---|---|
| six bedpe readers | `read_tempo_sv(detail =)` |
| six non-silent MAF filters | `non_silent()` |
| five FACETS amp/del rules | `cn_call(preset =)` |
| four `get_sample_data()` copies | `read_sample_data()` + `tempo_cohort_files()` |
| five chromosome conventions | `normalize_chrom()`, `chrom_factor()`, `autosomes()` |
| PCT denominators per report | `recurrence_table(denominator =)`, recorded in the DataDictionary |
| `cc()`, `len()`, `DATE()`, `write_xlsx()`, `halt()` from `.Rprofile` | `str_c()`, `length()`, `date_stamp()`, `write_report_xlsx()`, `wca_abort()` |
