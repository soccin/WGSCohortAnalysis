# CLAUDE.md

Notes for working on this repository. `README.md` is the reference,
`QUICKSTART.md` is the path from nothing to a cohort report.

## Patient data never enters git

The manifests name samples, and every workbook and figure derived from them
carries those names. None of it goes into a commit, a commit message, an
issue, or anything sent off the machine. Project numbers are de-identified
and are fine. Both shipped `.gitignore` files enforce this: the toolkit
ignores manifests and `data/raw/`, the project template ignores all of
`data/` except its README and all of `results/`.

Keep the manifests beside a clone, never inside one.

## What this is

A sourced R toolkit, not an R package. There is no `DESCRIPTION`, so
`testthat::test_local()` and `devtools::*` do not work here.

```
Rscript tests/run_tests.R          # unit tests + an end-to-end run on fixtures
WCA_SKIP_E2E=1 Rscript tests/run_tests.R   # unit tests only, much faster
```

Load it with `source(file.path(Sys.getenv("WCA_HOME"), "load.R"))`, which
attaches the dependencies, sources `R/*.R` in name order and sets
`wca.home`. Nothing in `R/` has top-level side effects.

```
load.R        wca_load()
R/            00 utils, 01-07 manifests/paths/readers (07 featureCounts
              expression), 10 cache, 20-21 genome and cn_call, 22-24
              events/recurrence/burden, 25 gene model and fusion exons,
              26 junction sequence triage, 27 alignment evidence (reads,
              soft clips, depth, insertion genotypes), 30-31 plots and
              BED, 40-41 xlsx writer and dictionary, 50 params and run axis
scripts/      the six pipeline stages, run against a project
bin/          wcaNewProject.R, wcaUpdateManifests.R
templates/    project scaffold copied by wcaNewProject.R
tests/        run_tests.R, testthat/, fixtures/miniCohort (renamed cell
              lines) and fixtures/miniModel (a hand-built GTF). The BAM
              and fasta fixtures for R/27 are synthesised in the test
              itself, so no sequencing data lives in the repository.
docs/         DATA_MODEL, METHODS, TEMPO_OUTPUT, CONSOLIDATION,
              CUSTOM_ANALYSIS, QUARTO_LINUX_INSTALL
ISSUES.md     known problems and open questions, one section each; read it
              before proposing work and delete a section when it is settled
```

An analysis is one folder: a fresh clone of this toolkit with the project
directory beside it, so nothing points outside that folder and the clone
records which toolkit produced the results.

## Rules that are easy to get wrong

- **CNV is autosomes-only.** Every CNV event, summary and plot drops X and Y
  through `filter_autosomes()`, because FACETS X calls are not sex-aware.
  Do not "fix" this. See `docs/METHODS.md`.
- **Sample identity is the tumor id.** `TID` in the manifests, `Sample` in
  the cohort table, with `NID` alongside. Manifest merges key on `TID`.
- **There are two manifests** because SNV and SV output can come from
  different Tempo runs. Expect them to disagree; the update report flags
  tumors that changed in only one.
- **Denominators are explicit.** `n_samples` for SNV and SV,
  `n_cnv_samples` (FACETS QC pass) for CNV.
- **Duplicate manifest rows.** `apply_cohort_rules()` collapses rows with
  the same md5 when the cohort is built. `--dedup` cleans the manifest
  itself and keeps the copy under `out/`, where Tempo writes. Rows with the
  same tumor but *different* md5 are reported, never resolved: that is a
  person's decision.
- **Chromosomes are unprefixed**, with FACETS 23/24 mapped to X/Y.
- **Never let a comma be read as a grouping mark.** `"90,0"` is two numbers,
  not 900, and no field in any file this toolkit reads ever means 900. But
  `readr::locale()` defaults `grouping_mark = ","`, so `read_tsv()`,
  `read_csv()`, `vroom()` and `type_convert()` all do exactly that, with no
  warning, and only when every value in the column happens to parse, so the
  same code is right on one cohort and wrong on the next. VCF INFO and FORMAT
  fields are full of comma-separated pairs: manta `PR`/`SR`, delly
  `CIPOS`/`CIEND`, per-allele `AC` and `AF`. `wca_load()` installs
  `wca_locale()` as the readr and vroom session default so this cannot
  happen, and `type_convert_silent()` passes it explicitly. Do not write a
  reader that guesses types under any other locale, do not pass
  `locale = readr::locale()` anywhere, and do not "fix" a column that came
  back as text by calling `parse_number()` on it. `tests/testthat/test-comma-numbers.R`
  guards this.
- **Never let YAML turn a letter into a boolean.** The yaml package follows
  YAML 1.1, where `y`, `n`, `yes`, `no`, `on` and `off` are booleans, so
  `ref_aa: Y` (tyrosine) reads as `TRUE`, `ref_aa: N` (asparagine) as
  `FALSE`, and `chrom: Y` as `TRUE`. `wca_read_params()` reads through
  `wca_read_yaml()`, where only `true` and `false` are logical, and stops
  if a key whose default is logical is written any other way. Read project
  yaml with these, never with `yaml::read_yaml()`; the one exception is the
  line in each stage script that finds the toolkit before it is loaded.
  `yaml::write_yaml()` writes `TRUE` as `yes`: pass
  `handlers = list(logical = yaml::verbatim_logical)`.
  `tests/testthat/test-params-yaml.R` guards this.

## Conventions

Tidyverse first: dplyr verbs, purrr over the apply family, the base pipe,
stringr and glue for strings, fs for paths, readr with
`show_col_types = FALSE`. Roxygen on every function. No trailing
whitespace.

Match the surrounding code. Add tests with a feature, and run the suite
before proposing a commit.

## Git

Work on a short-lived branch (`feat/manifest-dedup-01`), merge to master.
Commit only when asked. Draft the message in `/tmp` for review first: 50
character subject, 72 character body, no emoji.
