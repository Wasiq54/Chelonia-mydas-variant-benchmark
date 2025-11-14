

#!/bin/bash
cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/ || exit 1

out="normal_snp_jaccard.tsv"
echo -e "ToolA\tToolB\tJaccard" > "$out"

for A in BCF_normal_SNPs.bed Deepvariant_normal_SNPs.bed Freebayes_normal_SNPs.bed Gatk_normal_SNPs.bed Varscan_normal_SNPs.bed
do
  for B in BCF_normal_SNPs.bed Deepvariant_normal_SNPs.bed Freebayes_normal_SNPs.bed Gatk_normal_SNPs.bed Varscan_normal_SNPs.bed
  do
    j=$(bedtools jaccard -a "$A" -b "$B" | awk 'NR==2{print $3}')
    echo -e "${A%.bed}\t${B%.bed}\t${j}" >> "$out"
  done
done

echo "Done: $out"
