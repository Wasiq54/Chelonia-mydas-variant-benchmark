# Methods (Tools, Versions, Commands)

## Tools & versions (fill exact numbers you used)
- DeepVariant v1.x.x
- GATK v4.x.x
- FreeBayes v1.3.x
- VarScan v2.4.x
- BCFtools v1.xx
- BWA-MEM v0.7.17; Bowtie2 v2.5.x; Minimap2 v2.2x
- RTG Tools v3.12.1; FastQC v0.11.9; MultiQC v1.15; fastp v0.23.x; seqkit v2.x

## QC
fastqc --threads 16 W1-A_1.fastq.gz W1-A_2.fastq.gz -o CH-NORMS1/FastqcReport
fastqc --threads 16 S1_1.fastq.gz S1_2.fastq.gz -o Ab-NormS2/FastqcReport
multiqc "/media/work/New Volume1/DataofDNA/" -o "/media/work/New Volume1/DataofDNA/summary/"

## Alignment (examples)
bwa mem -t 32 rCheMyd1.pri.v2.fasta W1-A_1.fastq.gz W1-A_2.fastq.gz | samtools sort -o W1-A.bam
samtools index W1-A.bam

## Variant calling (examples)
# BCFtools
bcftools mpileup -Ou -f rCheMyd1.pri.v2.fasta W1-A.bam | bcftools call -mv -Oz -o W1-A.vcf.gz
tabix -p vcf W1-A.vcf.gz

# DeepVariant (sketch; replace with your exact command)
run_deepvariant --model_type=WGS --ref rCheMyd1.pri.v2.fasta --reads W1-A.bam --output_vcf W1-A.dv.vcf.gz --num_shards 32

## Filtering (if applied)
bcftools filter -e 'QUAL<30 || DP<10' -Oz -o W1-A.filtered.vcf.gz W1-A.vcf.gz
tabix -p vcf W1-A.filtered.vcf.gz

## Benchmark (RTG vcfeval example)
rtg vcfeval -b Truth.vcf.gz -c W1-A.filtered.vcf.gz -t ref.sdf -o rtg_out --vcf-score-field QUAL
