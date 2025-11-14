#!/bin/bash

# Input folder
INPUT="/after_filter"

# Output folder
OUT_SNP="/home/snps"

mkdir -p "$OUT_SNP"

cd "$INPUT"

for f in *.filtered.vcf.gz
do
    base=$(basename "$f" .filtered.vcf.gz)

    echo "Processing SNPs for: $f"

    bcftools view -v snps "$f" -Oz -o "$OUT_SNP/${base}_SNPs.vcf.gz"
    bcftools index "$OUT_SNP/${base}_SNPs.vcf.gz"
done

echo "✔ SNP files created in: $OUT_SNP"
