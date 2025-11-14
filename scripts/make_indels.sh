#!/bin/bash

# Input folder
INPUT="/home/after_filter"

# Output folder
OUT_INDEL="/home/indel"

mkdir -p "$OUT_INDEL"

cd "$INPUT"

for f in *.filtered.vcf.gz
do
    base=$(basename "$f" .filtered.vcf.gz)

    echo "Processing INDELs for: $f"

    bcftools view -v indels "$f" -Oz -o "$OUT_INDEL/${base}_INDELs.vcf.gz"
    bcftools index "$OUT_INDEL/${base}_INDELs.vcf.gz"
done

echo "✔ INDEL files created in: $OUT_INDEL"
