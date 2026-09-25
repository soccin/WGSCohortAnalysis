# QUICKSTART

One analysis lives in one folder: a fresh clone of the toolkit, with the
project directory beside it. Full detail is in `README.md`.

## 0. Clone the toolkit

```bash
cd ~/Work/MyAnalysis          # wherever this analysis is going to live
git clone git@github.com:soccin/WGSCohortAnalysis.git
export WCA_HOME=$PWD/WGSCohortAnalysis
```

Every command below uses `$WCA_HOME`, so nothing points outside this folder.
The clone is the record of which toolkit produced the results, so leave it
alone once a run starts and clone again for the next analysis.

## 1. Scaffold the project

```bash
Rscript $WCA_HOME/bin/wcaNewProject.R MyCohort myCohort
cd MyCohort
```

The directory must not exist. The second argument is the prefix on every
deliverable, as in `myCohort_CohortEvents_260910.xlsx`; it defaults to the
directory name. The scaffolder writes the clone's path into `toolkit:` in
`00.PARAMS.yml`, so the project always runs against the copy it was made
with.

## 2. Put the manifests in place

The manifests carry PHI. They name samples and the Tempo files behind them,
so they never go to GitHub: not in the clone, not in a project repository,
not in an issue. Keep them beside the clone and copy them in. The same goes
for everything under `results/`, which carries the same identifiers. Both
`.gitignore` files that ship with the toolkit already exclude them.

```bash
cp $SOMEWHERE/tempoSNVManifest_<stamp>.csv $SOMEWHERE/tempoSVManifest_<stamp>.csv data/raw/
```

Starting a cohort with no manifest at all? Build one by scanning the Tempo
output roots:

```bash
Rscript -e 'source(file.path(Sys.getenv("WCA_HOME"), "load.R")); \
  manifest_from_scan(scan_tempo_outputs("/path/to/Proj_NNNNN/out"), "data/raw")'
```

## 3. Edit 00.PARAMS.yml

Only these usually change:

| Key | What to set |
|---|---|
| `manifests.snv` / `manifests.sv` | the filenames you just copied |
| `cohort.include_projects` | the ProjNo values for this cohort, e.g. `[17495_I]`; leave `projects_file: null` (to select by Ucode instead, see `docs/METHODS.md`, Cohort selection) |
| `cohort.exclude_tid_regex` | drop cell lines and the like, e.g. `"_CL"` |
| `cohort.require_both` | `false` if a tumor is missing from one manifest |
| `report.genes_of_interest` | exact symbols, e.g. `[TP53, IKZF2]` |
| `run` | rename for a second pass with different cutoffs |

Leave `cnv.autosomes_only: true`. FACETS X and Y calls are not sex-aware.

## 4. Run

```bash
Rscript run_all.R            # all six stages
Rscript run_all.R 03 04      # rerun single stages after a parameter change
```

Stages are cohort, read, events, summaries, figures, report. Parsed Tempo
files are cached in `cache/reads/` and shared by every run, so a second run
is fast. Deliverables land in `results/<run>/`: two xlsx workbooks, a figures
pdf with PNGs, `tables/cohort.xlsx`, and `TOOLKIT_VERSION` plus a copy of the
parameters for provenance.

Check `results/<run>/tables/cohort.xlsx` first: it tells you which tumors made
the cohort and whether FACETS was found for each.

## 5. When new Tempo runs finish

List the output roots one per line, `#` comments allowed, optional `,ProjNo`
to override the project number:

```bash
cat > newProjects <<'EOF'
/path/to/Proj_18645_C/out
EOF

Rscript $WCA_HOME/bin/wcaUpdateManifests.R \
  data/raw/tempoSNVManifest_<stamp>.csv data/raw/tempoSVManifest_<stamp>.csv \
  newProjects data/raw --dedup
```

New tumors are appended, reruns of known tumors are updated in place, and
everything else is kept. Read `manifestUpdate_<stamp>.md` and the matching
csv before pointing `00.PARAMS.yml` at the new pair: they list every tumor
added, updated or dropped, with old and new values side by side. `--dedup`
drops rows naming a file already listed for the same tumor, keeping the copy
under `out/`.

## Gotchas

- The updater will not overwrite existing outputs. Use a new `--stamp` or
  `--force`.
- A tumor listed with two *different* files is reported, never resolved. Pick
  one by hand.
- `require_both: true` errors when the SNV and SV tumor sets differ. The
  report names the offenders.
- Keep manifests and their `attic` outside the clone; deleting a clone must
  never take data with it.
- Tests, after changing the toolkit: `Rscript $WCA_HOME/tests/run_tests.R`.
