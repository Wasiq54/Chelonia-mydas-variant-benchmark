

#!/bin/bash

# Folder containing INDEL BED files
cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/ || exit 1

# Abnormal INDEL BED file basenames (without .bed)
tools=(
  BCF_abnormal_INDELs
  Deepvariant_abnormal_INDELs
  Freebayes_abnormal_INDELs
  Gatk_abnormal_INDELs
  Varscan_abnormal_INDELs
)

# Output table
out="abnormal_indel_overlap.tsv"
echo -e "ToolA\tToolB\tOverlap" > "$out"

# Pairwise intersections
for ((i=0; i<${#tools[@]}; i++)); do
  for ((j=i+1; j<${#tools[@]}; j++)); do
    A="${tools[$i]}.bed"
    B="${tools[$j]}.bed"
    overlap=$(bedtools intersect -u -a "$A" -b "$B" | wc -l)
    echo -e "${tools[$i]}\t${tools[$j]}\t${overlap}" >> "$out"
  done
done

echo "Done. Results saved to $out"
