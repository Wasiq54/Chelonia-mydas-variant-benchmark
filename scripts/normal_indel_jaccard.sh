#!/bin/bash
cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/ || exit 1

out="normal_indel_jaccard.tsv"
echo -e "ToolA\tToolB\tJaccard" > "$out"

for A in BCF_normal_INDELs.bed Deepvariant_normal_INDELs.bed Freebayes_normal_INDELs.bed Gatk_normal_INDELs.bed Varscan_normal_INDELs.bed
do
  for B in BCF_normal_INDELs.bed Deepvariant_normal_INDELs.bed Freebayes_normal_INDELs.bed Gatk_normal_INDELs.bed Varscan_normal_INDELs.bed
  do
    j=$(bedtools jaccard -a "$A" -b "$B" | awk 'NR==2{print $3}')
    echo -e "${A%.bed}\t${B%.bed}\t${j}" >> "$out"
  done
done

echo "Done: $out"
