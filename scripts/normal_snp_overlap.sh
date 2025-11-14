


cd /home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/

tools=(BCF_normal_SNPs Deepvariant_normal_SNPs Freebayes_normal_SNPs Gatk_normal_SNPs Varscan_normal_SNPs)

echo -e "ToolA\tToolB\tOverlap" > normal_snp_overlap.tsv

for ((i=0; i<${#tools[@]}; i++)); do
    for ((j=i+1; j<${#tools[@]}; j++)); do
        A=${tools[$i]}.bed
        B=${tools[$j]}.bed
        overlap=$(bedtools intersect -u -a $A -b $B | wc -l)
        echo -e "${tools[$i]}\t${tools[$j]}\t$overlap" >> normal_snp_overlap.tsv
    done
done
