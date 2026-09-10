# Bundled reference data

All files here are read by `genome_info()` and the readers through
`wca_file("data", ...)`. Nothing else in the toolkit hard-codes a path.

| File | Purpose | Provenance |
|---|---|---|
| `cytoBand_hg19.txt.gz` | UCSC hg19 cytoBand table (chrom, start, end, band, stain). Chromosome lengths, centromeres (`acen` bands), arm boundaries and the acrocentric list are all derived from it, so it is the only genome file needed. | UCSC Table Browser, `hg19` / `cytoBand`, `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz`. Copied from `Analysis/260719/data/cytoBand_hg19.txt.gz`. |
| `sv_report_cols.txt` | Priority column order for the SV table (`read_tempo_sv()` puts these first). | `Analysis/AllTCells_260416/R/reportCols01`. |
| `maf_standard_cols.txt` | The MAF columns kept by `read_tempo_maf()` by default. Selected by name so both the 259 and 270 column Tempo schemas work. | Chosen from the Tempo `somatic.final.maf` header; the 11 `neo_*` columns are the only difference between schemas and are not in the list. |

## Coordinates

Tempo output for these projects is b37 (MAF `NCBI_Build == "GRCh37"`,
unprefixed chromosome names, FACETS codes X as 23). The cytoBand file uses
`chr` prefixes and is normalized on load by `normalize_chrom()`.

## Adding hg38

1. Download the hg38 cytoBand table to `data/cytoBand_hg38.txt.gz`.
2. Add an entry to the `genomes` list at the top of `R/20_genome.R` with the
   file name and the `NCBI_Build` aliases (`GRCh38`, `hg38`).
3. Every annotator already takes `genome =`, so no other code changes are
   needed. Set `genome: hg38` in a project `00.PARAMS.yml`.
