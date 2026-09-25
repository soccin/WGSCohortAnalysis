# ISSUES

Known problems and open questions in this toolkit, one section each. Not a
tracker: a short list of things a person has to decide or fix, with enough
context that the next reader does not have to rediscover them. Delete a
section when it is resolved and say so in the commit message.

## Cohort selection fails open without a warning

`scripts/01_cohort.R` unions `include_projects` with the workbook rows
`include_ucode` selects, and applies no project filter when the union is
empty. Three settings therefore select far more than intended and run to
completion (table in `docs/METHODS.md`, Cohort selection):

- `projects_file` with `include_ucode: []` selects every workbook project,
  whatever `include_projects` says. The APTL setup preview on 24 September
  2026 had this: 21 projects and 212 tumors instead of 5 and 22.
- `include_ucode` without `projects_file` is ignored.
- A `Ucode` that matches no workbook row selects nothing, so the union is
  empty and every project in the manifests is kept.

The docs now say to use one of the two ways. A guard in stage 01 would
make it impossible to get wrong: stop when `projects_file` is set and
`include_ucode` is empty, when `include_ucode` is set without
`projects_file`, and when `include_ucode` matches no workbook row; and
print the selected `ProjNo` values, not only their count.

Undecided. Raised 25 September 2026.

## delly_CIPOS and delly_CIEND are carried but never used

`read_tempo_sv()` pivots every delly FORMAT field, so the confidence
interval around each breakend arrives as `delly_CIPOS` and `delly_CIEND`.
Nothing in `R/` or `scripts/` reads either one, and neither reaches a
workbook sheet. They are two more columns in every SV table for no benefit.

They are comma-separated pairs, `-50,50`, which is how they came up: they
were two of the four columns the grouping-mark parsing bug corrupted. That
bug is fixed and they now parse correctly, as text, so this is not a
correctness question. It is only whether the columns are worth their width.

Three options. Drop them in the pivot. Keep them and split each into two
integer columns, which is what a consumer would want anyway. Leave them
alone on the grounds that a raw passthrough column costs little and
someone may want the interval later.

Undecided. Raised 12 September 2026.
