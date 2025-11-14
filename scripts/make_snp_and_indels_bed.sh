 mkdir -p snps_bed

for f in snps/*SNPs.vcf.gz
do
    base=$(basename $f .vcf.gz)
    bcftools query -f '%CHROM\t%POS\t%END\n' $f > snps_bed/${base}.bed
done

 


mkdir -p indels_bed

for f in indel/*INDELs.vcf.gz
do
    base=$(basename $f .vcf.gz)
    bcftools query -f '%CHROM\t%POS\t%END\n' $f > indels_bed/${base}.bed
done
