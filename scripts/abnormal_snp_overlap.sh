

#!/bin/bash

# Folder containing SNP BED files
cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/ || exit 1

# Abnormal SNP BED file basenames (without .bed)
tools=(
  BCF_abnormal_SNPs
  Deepvariant_abnormal_SNPs
  Freebayes_abnormal_SNPs
  Gatk_abnormal_SNPs
  Varscan_abnormal_SNPs
)

# Output table
out="abnormal_snp_overlap.tsv"
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
