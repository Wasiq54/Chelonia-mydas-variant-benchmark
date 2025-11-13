#!/bin/bash

OUTDIR="/home/All_tools_filtered_Varaints"
mkdir -p "$OUTDIR"

########################################
# DeepVariant
# AD = FORMAT (ref,alt), DP = FORMAT
########################################
for f in Deepvariant_*.vcf.gz; do
    [ -e "$f" ] || continue
    echo "Processing DeepVariant file: $f"

    base=$(basename "$f" .vcf.gz)
    outfile="$OUTDIR/${base}.filtered.vcf.gz"

    bcftools view \
      -i '(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / FMT/DP[0] > 0.02' \
      "$f" -Oz -o "$outfile"

    bcftools index "$outfile"

    count=$(bcftools view "$outfile" -H | wc -l)
    echo " → Output: $outfile"
    echo " → Variants passing filter: $count"
    echo
done


########################################
# BCFtools
# AD = FORMAT (ref,alt), DP = INFO
########################################
for f in BCF_*.vcf.gz; do
    [ -e "$f" ] || continue
    echo "Processing BCFtools file: $f"

    base=$(basename "$f" .vcf.gz)
    outfile="$OUTDIR/${base}.filtered.vcf.gz"

    bcftools view \
      -i '(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / INFO/DP > 0.02' \
      "$f" -Oz -o "$outfile"

    bcftools index "$outfile"

    count=$(bcftools view "$outfile" -H | wc -l)
    echo " → Output: $outfile"
    echo " → Variants passing filter: $count"
    echo
done


########################################
# FreeBayes
# AD = FORMAT (ref,alt), DP = FORMAT
########################################
for f in Freebayes_*.vcf.gz; do
    [ -e "$f" ] || continue
    echo "Processing FreeBayes file: $f"

    base=$(basename "$f" .vcf.gz)
    outfile="$OUTDIR/${base}.filtered.vcf.gz"

    bcftools view \
      -i '(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / FMT/DP[0] > 0.02' \
      "$f" -Oz -o "$outfile"

    bcftools index "$outfile"

    count=$(bcftools view "$outfile" -H | wc -l)
    echo " → Output: $outfile"
    echo " → Variants passing filter: $count"
    echo
done


########################################
# GATK (GenotypeGVCFs output)
# AD = FORMAT (ref,alt), DP = FORMAT
########################################
for f in Gatk_*.vcf.gz; do
    [ -e "$f" ] || continue
    echo "Processing GATK file: $f"

    base=$(basename "$f" .vcf.gz)
    outfile="$OUTDIR/${base}.filtered.vcf.gz"

    bcftools view \
      -i '(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / FMT/DP[0] > 0.02' \
      "$f" -Oz -o "$outfile"

    bcftools index "$outfile"

    count=$(bcftools view "$outfile" -H | wc -l)
    echo " → Output: $outfile"
    echo " → Variants passing filter: $count"
    echo
done


########################################
# VarScan
# AD = FORMAT (variant depth, scalar), DP = FORMAT
########################################
for f in Varscan_*.vcf.gz Varscan_*.reheader.vcf.gz; do
    [ -e "$f" ] || continue
    echo "Processing VarScan file: $f"

    base=$(basename "$f" .vcf.gz)
    base=${base%.reheader}
    outfile="$OUTDIR/${base}.filtered.vcf.gz"

    bcftools view \
      -i '(FILTER="PASS" || FILTER=".") && FMT/AD[0] > 3 && FMT/AD[0] / FMT/DP[0] > 0.02' \
      "$f" -Oz -o "$outfile"

    bcftools index "$outfile"

    count=$(bcftools view "$outfile" -H | wc -l)
    echo " → Output: $outfile"
    echo " → Variants passing filter: $count"
    echo
done
