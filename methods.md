<img width="591" height="170" alt="image" src="https://github.com/user-attachments/assets/1b8cb5ae-d194-4b46-9904-fe9242839a8e" /># Methods (Tools, Versions, Commands)

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



**##Trimming**
Five widely used read-preprocessing tools were evaluated for quality and adapter trimming of Chelonia mydas paired-end Illumina NovaSeq reads (151 bp). All runs were performed using 20 threads on a high-performance workstation (Intel Xeon, 512 GB RAM, Ubuntu 20.04 LTS).

Tool Versions and Sources:
------------------------------------------------------------
1. fastp v0.26.0 — OpenGene GitHub
   Fast, all-in-one trimmer with automatic adapter detection, base correction, and per-sample JSON/HTML QC reports (released August 2024).

2. Trimmomatic v0.39 — Trimmomatic GitHub
   Standard Java-based trimmer using sliding-window and adapter clipping; stable release since 2019.

3. BBDuk (BBMap suite) v39.01 — BBMap SourceForge
   Adapter and quality trimming based on k-mer matching; included in BBTools v39.01 (February 2024).

4. Trim Galore! v0.6.8 — Babraham Bioinformatics
   Wrapper around Cutadapt; performs quality trimming and adapter removal (December 2023).

5. Cutadapt v5.2 — Bioconda / GitHub
   Python-based adapter trimmer supporting multithreading (released 2025; Python 3.12 compatible).
------------------------------------------------------------

Workflow and Rationale
All five trimmers were applied to both the normal (CH-NORMS1) and abnormal (Ab-NormS2) datasets to compare performance in adapter removal, read retention, and GC content preservation.
Among them, fastp (v0.26.0) demonstrated the best overall balance of speed, accuracy, and comprehensive reporting. Therefore, fastp-trimmed reads were selected for all downstream analyses, including alignment, variant calling, and benchmarking.

The exact shell commands and parameters for each tool are provided in:
scripts/commands_trimming_all_tools.txt

Quality improvements were confirmed using FastQC (v0.11.9) and summarized with MultiQC (v1.15) after trimming.





## Alignment (examples)

Alignment
------------------------------------------------------------
Three mappers were evaluated for read alignment against the Chelonia mydas reference genome (rCheMyd1.pri.v2):
1. BWA-MEM v0.7.17
2. Bowtie2 v2.5.4
3. Minimap2 v2.28

Each aligner was executed using 15–40 threads on paired-end 151 bp Illumina NovaSeq reads.
Outputs were sorted and indexed with Samtools (v1.22.1).
Exact shell commands are provided in scripts/commands_alignment.txt.
Mapping quality was evaluated using Samtools flagstat and alignment summary metrics.







Variant Calling
------------------------------------------------------------
Five variant calling tools were used to detect SNPs and indels from BWA-MEM–aligned reads:
DeepVariant v1.9.0, GATK v4.6.2.0, FreeBayes v1.3.10, VarScan v2.4.6, and BCFtools v1.22.
Each tool was executed on both normal and abnormal Chelonia mydas datasets using the same reference genome (rCheMyd1.pri.v2).
Command-line parameters for reproducibility are provided in scripts/commands_variant_calling.txt.




















## Filtering (if applied)
bcftools filter -e 'QUAL<30 || DP<10' -Oz -o W1-A.filtered.vcf.gz W1-A.vcf.gz
tabix -p vcf W1-A.filtered.vcf.gz

## Benchmark (RTG vcfeval example)
rtg vcfeval -b Truth.vcf.gz -c W1-A.filtered.vcf.gz -t ref.sdf -o rtg_out --vcf-score-field QUAL
