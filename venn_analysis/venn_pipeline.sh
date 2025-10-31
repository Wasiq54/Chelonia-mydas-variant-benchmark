#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Five-caller SNP overlap prep for Venn diagrams (Normal/Abnormal)
# Requires: bcftools (>=1.10), tabix
# Author: Wasiq Aslam
# ============================================================

# --- EDIT THESE PATHS IF NEEDED --------------------------------------------
REFERENCE="/media/fahad/DNA_Work2/validation2/ref/reference_genome.fasta"

# Normal and Abnormal working dirs (where outputs will be written)
NORMAL_DIR="/media/fahad/DNA_Work2/concordance/ForvanDiagram"
ABNORM_DIR="/media/fahad/DNA_Work2/concordance/ForvanDiagram/Abnormal"
# ---------------------------------------------------------------------------

# Ensure dirs exist
mkdir -p "$NORMAL_DIR" "$ABNORM_DIR"

# ----------------------------- NORMAL --------------------------------------
echo "==> NORMAL dataset"
cd "$NORMAL_DIR"

declare -A NORMAL_VCFS=(
  ["DeepVariant"]="../DeepVariant_NORM_Sorted_output.vcf.gz"
  ["GATK"]="../GatK_CH-NORMS1_raw_variants.vcf.gz"
  ["FreeBayes"]="../freeBayes_CH-NORMS1_variants.vcf.gz"
  ["VarScan"]="../Varscan_Normal.vcf.gz"
  ["BCFTools"]="../BCFTools_CH-NORMS1_variants.vcf.gz"
)

for tool in DeepVariant GATK FreeBayes VarScan BCFTools; do
  in_vcf="${NORMAL_VCFS[$tool]}"
  echo "  • $tool : $in_vcf"

  bcftools norm -f "$REFERENCE" -m -both "$in_vcf" -Oz -o "${tool}_Normal_normalized.vcf.gz"
  bcftools sort "${tool}_Normal_normalized.vcf.gz" -Oz -o "${tool}_Normal_sorted.vcf.gz"
  bcftools index -f "${tool}_Normal_sorted.vcf.gz"

  # SNP-only coordinate list for Venn
  bcftools view -v snps "${tool}_Normal_sorted.vcf.gz" \
    | bcftools query -f '%CHROM:%POS\n' \
    | sort -u > "${tool}_Normal_SNPs.txt"
done

# ---------------------------- ABNORMAL -------------------------------------
echo "==> ABNORMAL dataset"
cd "$ABNORM_DIR"

declare -A ABNORM_VCFS=(
  ["DeepVariant"]="../../../DeepVariant_Abnorm_sorted_output.vcf.gz"
  ["GATK"]="../../../GatK_Ab-NormS2_raw_variants.vcf.gz"
  ["FreeBayes"]="../../../freeBayes_Ab-NormS2_variants.vcf.gz"
  ["VarScan"]="../../../Varscan_Ab-NormS2.vcf.gz"
  ["BCFTools"]="../../../BCFTools_Ab-NormS2_variants.vcf.gz"
)

for tool in DeepVariant GATK FreeBayes VarScan BCFTools; do
  in_vcf="${ABNORM_VCFS[$tool]}"
  echo "  • $tool : $in_vcf"

  bcftools norm -f "$REFERENCE" -m -both "$in_vcf" -Oz -o "${tool}_Abnorm_normalized.vcf.gz"
  bcftools sort "${tool}_Abnorm_normalized.vcf.gz" -Oz -o "${tool}_Abnorm_sorted.vcf.gz"
  bcftools index -f "${tool}_Abnorm_sorted.vcf.gz"

  bcftools view -v snps "${tool}_Abnorm_sorted.vcf.gz" \
    | bcftools query -f '%CHROM:%POS\n' \
    | sort -u > "${tool}_Abnorm_SNPs.txt"
done

echo "Done. SNP lists written to:"
echo "  $NORMAL_DIR/*_Normal_SNPs.txt"
echo "  $ABNORM_DIR/*_Abnorm_SNPs.txt"
