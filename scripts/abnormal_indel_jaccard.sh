#!/bin/bash
cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/ || exit 1

out="abnormal_indel_jaccard.tsv"
echo -e "ToolA\tToolB\tJaccard" > "$out"

for A in BCF_abnormal_INDELs.bed Deepvariant_abnormal_INDELs.bed Freebayes_abnormal_INDELs.bed Gatk_abnormal_INDELs.bed Varscan_abnormal_INDELs.bed
do
  for B in BCF_abnormal_INDELs.bed Deepvariant_abnormal_INDELs.bed Freebayes_abnormal_INDELs.bed Gatk_abnormal_INDELs.bed Varscan_abnormal_INDELs.bed
  do
    j=$(bedtools jaccard -a "$A" -b "$B" | awk 'NR==2{print $3}')
    echo -e "${A%.bed}\t${B%.bed}\t${j}" >> "$out"
  done
done

echo "Done: $out"
