# <project name>

WGS cohort analysis of Tempo output, built with
[WGSCohortAnalysis](../WGSCohortAnalysis). Layout follows
https://github.com/soccin/r-proj-org.

```
00.PARAMS.yml        parameters for this project (run name, manifests, cutoffs)
run_all.R            Rscript run_all.R runs stages 01..06 from the toolkit
data/raw/            the two Tempo manifests (+ optional projects xlsx)
cache/reads/         parsed Tempo files, content-addressed (shared by all runs)
cache/<run>/         per-stage intermediate rds
results/<run>/       workbooks, figures, TOOLKIT_VERSION, PARAMS copy
analysis/            report sources (.qmd, .Rmd); never in the root
```

## Run

```
export WCA_HOME=$PWD/../WGSCohortAnalysis    # the clone beside this project
Rscript run_all.R
```

Re-running is cheap: stage 02 hits `cache/reads/`. Change `run:` in
`00.PARAMS.yml` to keep a previous run's results.

## Reports

A Quarto report lives in `analysis/`, reads the stage cache, and renders
into the run's results folder. Build paths in it from `here::here()`;
the code runs in `analysis/`. Render from inside `analysis/` (from the
project root, Quarto 1.10 fails to find its support files):

```
cd analysis
quarto render reportNN_<topic>.qmd -P run:run01 \
  --output-dir ../results/run01 --output <prefix>_<Topic>_<yymmdd>.html
```

## Custom analyses

```r
source(file.path(Sys.getenv("WCA_HOME"), "load.R"))
PARAMS <- wca_read_params("00.PARAMS.yml")
snv <- stage_load(PARAMS, "03_events", "snv")
```

See `docs/CUSTOM_ANALYSIS.md` in the toolkit for worked examples.
