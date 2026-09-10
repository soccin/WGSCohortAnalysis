# Test fixtures

`miniCohort/` is a trimmed Tempo run with three tumor/normal pairs
(CL01__CL01N, CL02__CL02N, CL03__CL03N) built by `make_fixtures.sh` from the
Proj_17297_B cell-line pairs (IRIS run 17329). Ids are renamed in file
names and contents.

Per pair: MAF (about 100 rows chosen to include LoF classes, missense,
silent, intronic, TERT and TP53 rows and chromosome X), bedpe (all header
lines and about 30 rows including translocations, an in-frame fusion and
self-pairs), FACETS gene_level (head plus TP53, IKZF2, GATA3, BTK,
chromosome 23 and non-diploid rows), arm_level, facets_qc, cncf and seg
heads, meta_data sample_data, and a synthesized `multiqc/QC_Status.txt`
(this Tempo vintage does not write one). `cohort_level/default_cohort/`
holds the concatenated sample_data and a synthetic alignment_qc.txt on the
real header.

Regenerate with `bash tests/fixtures/make_fixtures.sh` (cluster paths).
