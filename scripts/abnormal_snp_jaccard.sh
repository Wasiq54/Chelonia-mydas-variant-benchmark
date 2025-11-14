#!/bin/bash
cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/ || exit 1

out="abnormal_snp_jaccard.tsv"
echo -e "ToolA\tToolB\tJaccard" > "$out"

for A in BCF_abnormal_SNPs.bed Deepvariant_abnormal_SNPs.bed Freebayes_abnormal_SNPs.bed Gatk_abnormal_SNPs.bed Varscan_abnormal_SNPs.bed
do
  for B in BCF_abnormal_SNPs.bed Deepvariant_abnormal_SNPs.bed Freebayes_abnormal_SNPs.bed Gatk_abnormal_SNPs.bed Varscan_abnormal_SNPs.bed
  do
    j=$(bedtools jaccard -a "$A" -b "$B" | awk 'NR==2{print $3}')
    echo -e "${A%.bed}\t${B%.bed}\t${j}" >> "$out"
  done
done

echo "Done: $out"
