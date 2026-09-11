# WGSCohortAnalysis

Cohort-level analysis of Tempo whole-genome output: SNV (MAF), CNV (FACETS)
and SV (bedpe) readers, annotators, recurrence summaries, figures and xlsx
reports, plus a standard six-stage pipeline that runs against a project's
`00.PARAMS.yml`. It replaces the per-project scripts surveyed in
`docs/CONSOLIDATION.md` with one set of functions and one set of rules
(`docs/METHODS.md`).

Version: see `VERSION`. Genome: hg19/b37 bundled; every annotator takes
`genome =` (see `data/README.md` to add hg38).

## Requirements

R >= 4.3 with tidyverse (dplyr, tidyr, purrr, readr, stringr, tibble,
forcats, ggplot2), vroom, fs, glue, rlang, yaml, readxl, openxlsx,
data.table, digest, patchwork, scales, testthat and withr (tests only).
No dependency on `~/.Rprofile`.

## Quick start: custom use

```r
source(file.path(Sys.getenv("WCA_HOME"), "load.R"))   # WCA_HOME = your clone of this toolkit

snv_m <- read_manifest("data/raw/tempoSNVManifest.csv") |> apply_cohort_rules(exclude_tid_regex = "_CL")
sv_m  <- read_manifest("data/raw/tempoSVManifest.csv")  |> apply_cohort_rules(exclude_tid_regex = "_CL")
cohort <- build_cohort(snv_m, sv_m)

maf <- read_cohort_maf(cohort, cache_dir = "cache/reads")       # non-silent rows, standard columns
sv  <- read_cohort_sv(cohort, cache_dir = "cache/reads") |> sv_events()
fac <- read_cohort_facets(cohort, cache_dir = "cache/reads")
qc  <- facets_qc_gate(fac$qc)

snv <- snv_events(maf)
cnv <- cnv_events(fac$gene, qc)
svg <- sv_gene_events(sv)

n <- nrow(cohort); n_cnv <- length(cnv_samples(qc))
tbl <- event_summary(list(
  recurrence_table(snv, n, "SNV"),
  recurrence_table(cnv, n_cnv, "CNV"),
  recurrence_table(svg, n, "SV")), n)
write_report_xlsx(list(AllEvents = tbl), "events.xlsx")
plot_oncoprint(bind_events(snv, cnv, svg))
```

`docs/CUSTOM_ANALYSIS.md` has four worked examples (a TP53 per-sample
report, a gene-set matrix, a cohort subset by regex, events shared by
sample groups).

## Quick start: standard run

One analysis is one folder: a clone of this toolkit with the project beside
it, so nothing points outside that folder and the clone records which
toolkit produced the results. `QUICKSTART.md` walks through it. The
manifests and everything under `results/` carry PHI and never go to a
remote; both shipped `.gitignore` files exclude them.

```
cd ~/Work/MyAnalysis
git clone git@github.com:soccin/WGSCohortAnalysis.git
export WCA_HOME=$PWD/WGSCohortAnalysis

Rscript $WCA_HOME/bin/wcaNewProject.R MyCohort myCohort
cd MyCohort
cp .../tempoSNVManifest_*.csv .../tempoSVManifest_*.csv data/raw/
$EDITOR 00.PARAMS.yml          # manifests, cohort rules, genes of interest, cutoffs
Rscript run_all.R               # or: Rscript run_all.R 03 04 to rerun stages
```

Outputs land in `results/<run>/`:

| File | Content |
|---|---|
| `<prefix>_CohortEvents_<yymmdd>.xlsx` | Samples, FacetsQC, AllEvents, LossOfFunc, SNVGenes, CNVGenes, CNVGains, CNVBands, CNVArms, SVGenes, SVPairs, Fusions, SVPartners, GeneSet, GeneSetLoF, DataDictionary |
| `<prefix>_GeneEvents_<yymmdd>.xlsx` | SNV, CNV, SV detail rows for the genes of interest |
| `figures/<prefix>_figures_<yymmdd>.pdf` + PNGs | oncoprint, recurrent genes, SV burden / mix / size, cytoband gain-loss, arm heatmap, purity-ploidy, TMB |
| `tables/cohort.xlsx` | the cohort table (and the projects sheet when used) |
| `TOOLKIT_VERSION`, `00.PARAMS.yml` | provenance: toolkit version + git SHA, parameters used |

Stage intermediates are rds files under `cache/<run>/NN_stage/`; parsed
Tempo files are cached by content under `cache/reads/` and shared by every
run (the one documented deviation from r-proj-org).

## Updating the manifests

When a project is rerun or new projects finish, list their Tempo output
roots (one per line, `#` comments allowed, optional `,ProjNo` to override
the project number derived from the run directory) and run

```
Rscript $WCA_HOME/bin/wcaUpdateManifests.R \
  tempoSNVManifest_260726.csv tempoSVManifest_260726.csv updatedOrNewProjects
```

It scans the roots for `somatic/<PAIR>/` MAF and bedpe files, merges them
into both manifests keyed on `TID` (new tumors appended, reruns of known
tumors updated in place, everything else kept) and writes
`tempoSNVManifest_<yymmdd>.csv`, `tempoSVManifest_<yymmdd>.csv`,
`manifestUpdate_<yymmdd>.md` and `manifestUpdate_<yymmdd>.csv` next to the
input manifests. The report lists every added and updated tumor with what
changed and flags tumors that changed in only one of the two manifests; the
csv is the same set of tumors as one flat table, one row per manifest and
`TID`, old and new values side by side, for review in a spreadsheet.

Tumors listed more than once are always reported, in the run log and in a
`duplicates` section of the report, split into two kinds: the same file
reached by two paths (one md5, safe to drop) and two different files for
one tumor (a decision about which run is current, which the toolkit will
not make for you). `--dedup` drops the first kind and every dropped row is
listed in the report and the csv. Of a redundant group the copy under
`out/` is kept, because that is where Tempo writes and a second path for
the same md5 is a copy of it; `--prefer=<regex on PATH>` changes that rule
and `--prefer=` disables it, leaving the choice to `--dedup` (first row of
the group), `--dedup=newest` or `--dedup=oldest` (by file age).
Note that `apply_cohort_rules()` already collapses same-md5 rows when the
cohort is built, so deduplicating the manifest is tidiness, not a fix for a
wrong analysis.

Options: `[out_dir]`, `--stamp=yymmdd`, `--force` (overwrite outputs),
`--dedup[=first|newest|oldest]`, `--prefer=<regex>`. The underlying functions are
`update_manifest()`, `manifest_update_changes()`, `manifest_duplicates()`
and `dedup_manifest()`.

## Layout

```
load.R          wca_load(): sources R/*.R, attaches deps, sets option wca.home
R/              00 utils, 01-06 manifests/paths/readers, 10 cache, 20-21 genome and
                cn_call, 22-24 events/recurrence/burden, 30-31 plots and BED,
                40-41 xlsx writer and dictionary, 50 params and run axis
scripts/        the six pipeline stages (run against a project)
templates/      project scaffold (00.PARAMS.yml, run_all.R, README, .gitignore)
bin/            wcaNewProject.R, wcaUpdateManifests.R
data/           cytoBand hg19, column lists (see data/README.md)
tests/          Rscript tests/run_tests.R ; fixtures/miniCohort (3 renamed cell-line pairs)
docs/           DATA_MODEL, METHODS, TEMPO_OUTPUT, CONSOLIDATION, CUSTOM_ANALYSIS,
                QUARTO_LINUX_INSTALL
```

## Conventions

Tidyverse style, base pipe, roxygen on every function, no top-level side
effects in `R/`. Sample identity is the tumor id (`Sample`) with `NID`
alongside. Chromosomes are unprefixed with FACETS 23/24 mapped to X/Y.
CNV analysis is autosomes-only by default: every CNV event, summary and
plot function drops X and Y through `filter_autosomes()` because FACETS
X calls are not sex-aware (see docs/METHODS.md). Denominators are
explicit: `n_samples` for SNV and SV, `n_cnv_samples` (FACETS QC pass)
for CNV.

## Tests

```
Rscript tests/run_tests.R
```

Runs the unit tests and an end-to-end run of `run_all.R` on the fixture
cohort (set `WCA_SKIP_E2E=1` to skip the latter).

## Not in this iteration

DAD1/BLAT read parsing, OGM/Delly/DRAGEN comparisons, CNV versus RNA,
IMPACT/cBioPortal extraction, invoice and clinical reformatters, coverage
QC triage (all listed with paths in `docs/CONSOLIDATION.md`); hg38 data,
ComplexHeatmap/circos, an HTML report (see `docs/QUARTO_LINUX_INSTALL.md`),
SV signatures, ClusterSV harvesting, mutual-exclusivity statistics. The
proposal list in `kej-AML_OGM/docs/PROPOSED-tempo-sv-analyses.md` is the
source for the next candidates.
