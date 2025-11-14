tools=(BCF_normal_INDELs Deepvariant_normal_INDELs Freebayes_normal_INDELs Gatk_normal_INDELs Varscan_normal_INDELs)

echo -e "ToolA\tToolB\tOverlap" > normal_indel_overlap.tsv

for ((i=0; i<${#tools[@]}; i++)); do
    for ((j=i+1; j<${#tools[@]}; j++)); do
        A=${tools[$i]}.bed
        B=${tools[$j]}.bed
        overlap=$(bedtools intersect -u -a $A -b $B | wc -l)
        echo -e "${tools[$i]}\t${tools[$j]}\t$overlap" >> normal_indel_overlap.tsv
    done
done
