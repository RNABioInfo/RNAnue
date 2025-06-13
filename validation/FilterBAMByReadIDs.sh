#!/bin/bash

set -euo pipefail

bam="$1"          # e.g. input.bam
tsv="$2"          # e.g. clusters.tsv
outdir="${3:-.}"  # optional output directory

mkdir -p "$outdir"

# skip header with tail, then for each line split on TAB into $label and $ids
tail -n +2 "$tsv" | \
while IFS=$'\t' read -r label ids; do
    [[ -z $label || -z $ids ]] && continue
    # turn commas into newlines on the fly and feed into -N
    samtools view -b -N <(tr ',' '\n' <<< "$ids") "$bam" \
    > "$outdir/${label}.bam"

    samtools sort -o "$outdir/${label}.bam" "$outdir/${label}.bam"

    samtools index "$outdir/${label}.bam"
done
