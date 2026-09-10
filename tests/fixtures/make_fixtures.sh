#!/usr/bin/env bash
# make_fixtures.sh
#
# Rebuild tests/fixtures/miniCohort from three real Tempo cell-line pairs
# (Proj_17297_B, IRIS run 17329). Files are trimmed to a few hundred rows and
# sample ids are renamed to CL01__CL01N, CL02__CL02N, CL03__CL03N. The
# QC_Status.txt files are synthesized because this Tempo vintage does not
# emit a per-pair multiqc directory. Run from anywhere:
#
#   SRC=<tempo_run_dir> OLD_N=<normal_sample_id> bash tests/fixtures/make_fixtures.sh
#
# The run directory and the normal sample id are passed in rather than
# recorded here: the run is only reachable on the BIC cluster and the id
# names a patient. The generated fixtures carry renamed ids only, and are
# committed, so the tests never need either.

set -euo pipefail

SRC=${SRC:?set SRC to the Tempo run directory to sample from (Proj_17297_B, IRIS run 17329)}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DST=$HERE/miniCohort/out/miniCohort
FAC=facets0.5.14c1000pc5000
OLD_N=${OLD_N:?set OLD_N to the normal sample id the three pairs share}

rm -rf "$HERE/miniCohort"
mkdir -p "$DST/somatic" "$DST/cohort_level/default_cohort"

declare -A PAIRS=(
  [p17297b_Myla]=CL01
  [p17297b_Seax]=CL02
  [p17297b_KEJ_H9]=CL03
)

# Column indices (1-based) in the Tempo MAF
#   1 Hugo_Symbol  5 Chromosome  9 Variant_Classification
trim_maf() {
  awk -F'\t' -v OFS='\t' '
    NR == 1 { print; next }
    /^#/ { next }
    $1 == "TERT" || $1 == "TP53" { print; next }
    $9 ~ /Frame_Shift|Nonsense|Splice_Site|In_Frame|Translation_Start|Nonstop/ { if (lof++ < 40) print; next }
    $5 == "X" { if (chrx++ < 4) print; next }
    $9 == "Missense_Mutation" { if (mis++ < 25) print; next }
    $9 == "Silent" { if (sil++ < 8) print; next }
    $9 == "Intron" { if (intr++ < 6) print; next }
    $9 == "IGR" { if (igr++ < 3) print; next }
    $9 == "5'"'"'Flank" { if (fl5++ < 3) print; next }
    $9 == "3'"'"'UTR" || $9 == "5'"'"'UTR" || $9 == "RNA" || $9 == "3'"'"'Flank" { if (oth++ < 4) print; next }
  ' "$1"
}

# Column indices in the bedpe: 11 TYPE, 27 gene1, 30 gene2, 33 fusion
trim_bedpe() {
  awk -F'\t' -v OFS='\t' '
    /^##/ { print; next }
    /^#/ { print; next }
    $33 ~ /Protein Fusion: in frame/ { if (fus++ < 3) print; next }
    $27 == $30 && $27 != "." { if (self++ < 3) print; next }
    $11 == "BND" { if (tra++ < 6) print; next }
    $11 == "DEL" { if (del++ < 7) print; next }
    $11 == "DUP" { if (dup++ < 5) print; next }
    $11 == "INV" { if (inv++ < 5) print; next }
    $11 == "INS" { if (ins++ < 2) print; next }
  ' "$1"
}

# Column indices in gene_level.txt: 2 gene, 3 chrom, 24 cn_state
trim_gene_level() {
  awk -F'\t' -v OFS='\t' '
    NR == 1 { print; next }
    $2 == "TP53" || $2 == "IKZF2" || $2 == "GATA3" || $2 == "BTK" { print; next }
    $3 == "23" { if (chrx++ < 8) print; next }
    $24 != "DIPLOID" { if (alt++ < 60) print; next }
    { if (dip++ < 60) print; next }
  ' "$1"
}

sample_data_header=""
for old in "${!PAIRS[@]}"; do
  new=${PAIRS[$old]}
  new_n="${new}N"
  old_pair="${old}__${OLD_N}"
  new_pair="${new}__${new_n}"
  ren() { sed -e "s/${old}/${new}/g" -e "s/${OLD_N}/${new_n}/g"; }

  S=$SRC/somatic/$old_pair
  D=$DST/somatic/$new_pair
  mkdir -p "$D/combined_mutations" "$D/combined_svs" "$D/facets/$new_pair/$FAC" \
           "$D/meta_data" "$D/multiqc"

  trim_maf "$S/combined_mutations/$old_pair.somatic.final.maf" | ren \
    > "$D/combined_mutations/$new_pair.somatic.final.maf"
  trim_bedpe "$S/combined_svs/$old_pair.final.bedpe" | ren \
    > "$D/combined_svs/$new_pair.final.bedpe"

  F=$S/facets/$old_pair
  ren < "$F/$old_pair.facets_qc.txt" > "$D/facets/$new_pair/$new_pair.facets_qc.txt"
  ren < "$F/${old_pair}_OUT.txt"     > "$D/facets/$new_pair/${new_pair}_OUT.txt"
  trim_gene_level "$F/$FAC/$old_pair.gene_level.txt" | ren \
    > "$D/facets/$new_pair/$FAC/$new_pair.gene_level.txt"
  ren < "$F/$FAC/$old_pair.arm_level.txt" > "$D/facets/$new_pair/$FAC/$new_pair.arm_level.txt"
  ren < "$F/$FAC/$old_pair.qc.txt"        > "$D/facets/$new_pair/$FAC/$new_pair.qc.txt"
  head -40 "$F/$FAC/${old_pair}_hisens.cncf.txt" | ren > "$D/facets/$new_pair/$FAC/${new_pair}_hisens.cncf.txt"
  head -40 "$F/$FAC/${old_pair}_hisens.seg"      | ren > "$D/facets/$new_pair/$FAC/${new_pair}_hisens.seg"
  head -40 "$F/$FAC/${old_pair}_purity.cncf.txt" | ren > "$D/facets/$new_pair/$FAC/${new_pair}_purity.cncf.txt"
  head -40 "$F/$FAC/${old_pair}_purity.seg"      | ren > "$D/facets/$new_pair/$FAC/${new_pair}_purity.seg"

  ren < "$S/meta_data/$old_pair.sample_data.txt" > "$D/meta_data/$new_pair.sample_data.txt"

  printf '\tStatus\tReason\n%s\tpass\tPassed Tempo Criteria\n%s\tpass\tPassed Tempo Criteria\n' \
    "$new" "$new_n" > "$D/multiqc/$new_pair.QC_Status.txt"

  if [ -z "$sample_data_header" ]; then
    head -1 "$D/meta_data/$new_pair.sample_data.txt" > "$DST/cohort_level/default_cohort/sample_data.txt"
    sample_data_header=done
  fi
  tail -n +2 "$D/meta_data/$new_pair.sample_data.txt" >> "$DST/cohort_level/default_cohort/sample_data.txt"
done

# A small alignment_qc.txt in the older cohort_level format (synthetic values
# on the real header) so read_alignment_qc() has something to parse.
{
  printf 'Sample\tTotalReads\tNumberDuplicateMarked\tFractionDuplicateMarked\tNumberMapped\tFractionMapped\tNumberUnmapped\tFractionUnmapped\tNumberSecondaryAlignments\tFractionSecondaryAlignments\tTotalPairs\tNumberMappedProperPairs\tFractionMappedProperPairs\tMedianReadLength\tMedianInsertSize\tMedianCoverage\n'
  for new in CL01 CL02 CL03; do
    printf '%s\t1500000000\t300000000\t0.2\t1200000000\t0.8\t300000\t0.0002\t0\t0\t750000000\t740000000\t0.98\t151:151\t400\t65\n' "$new"
    printf '%sN\t1200000000\t200000000\t0.17\t1000000000\t0.83\t160000\t0.0001\t0\t0\t600000000\t590000000\t0.98\t151:151\t399\t50\n' "$new"
  done
} > "$DST/cohort_level/default_cohort/alignment_qc.txt"

du -sh "$HERE/miniCohort"
find "$HERE/miniCohort" -type f | sort
