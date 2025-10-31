# Chelonia-mydas-variant-benchmark

**Benchmarking DeepVariant and Conventional Variant Callers for the Chelonia mydas Genome** — including scripts, workflows, and reproducibility materials.

This repository contains all command scripts and workflows used in:
"Benchmarking DeepVariant and Conventional Variant Callers for Consensus-Based Analysis of the Chelonia mydas Genome."

---

## Data Availability

The whole-genome sequencing datasets analyzed in this study are publicly available in the NCBI Sequence Read Archive (SRA) under BioProject accession:

**BioProject:** [PRJNA1312993](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1312993)

- **Normal sample:** BioSample [SAMN50881756](https://www.ncbi.nlm.nih.gov/biosample/SAMN50881756), SRA run accession [SRR35202080](https://www.ncbi.nlm.nih.gov/sra/SRR35202080)
- **Abnormal sample:** BioSample [SAMN50881755](https://www.ncbi.nlm.nih.gov/biosample/SAMN50881755), SRA run accession [SRR35202079](https://www.ncbi.nlm.nih.gov/sra/SRR35202079)

Processed data, including trimmed reads, alignment (BAM), and variant call (VCF) files, are available upon reasonable request from the corresponding author.

---

## Workflow Overview

1. **Quality Control:**
   - Tools: FastQC and MultiQC
   - Commands: See `scripts/commands_alignment.txt`

2. **Read Trimming:**
   - Tools evaluated: fastp, Trimmomatic, TrimGalore, Cutadapt, BBDuk
   - Main trimmer used: **fastp v0.26.0**
   - Commands: `scripts/commands_trimming_all_tools.txt`

3. **Alignment:**
   - Primary mapper: **BWA-MEM**
   - Comparative mappers: Bowtie2 and Minimap2
   - Commands: `scripts/commands_alignment.txt`

4. **Variant Calling:**
   - Tools: DeepVariant, GATK, FreeBayes, VarScan, BCFtools
   - Commands: `scripts/commands_variant_calling.txt`

5. **Consensus and Benchmarking:**
   - Intersection filtering: `bcftools isec` (3+, 4+, 5+ thresholds)
   - Benchmarking: **RTG vcfeval**
   - Commands: `scripts/commands_consensus_truthset.txt`

6. **Concordance and Discordance Visualization:**
   - Venn diagram generation: `venn_analysis/venn_pipeline.sh`
   - Plotting scripts: `venn_plot_normal.R` and `venn_plot_abnormal.R`
   - Missed/unique SNP detection: `discordance_analysis.py`

---

## Software Versions

All analyses were performed using the latest stable releases (as of 2025) in Conda or Docker environments.

- DeepVariant v1.9.0
- GATK v4.6.2.0
- FreeBayes v1.3.10
- VarScan v2.4.6
- BCFtools v1.22
- Samtools v1.22.1
- BWA-MEM v0.7.17
- Bowtie2 v2.5.2
- Minimap2 v2.28
- RTG Tools v3.12.1
- FastQC v0.11.9
- MultiQC v1.15
- fastp v0.26.0
- seqkit v2.8.2

---




## Reproducibility and Citation

All workflows are detailed in `methods.md` with complete parameters and environment details.  
Analyses were executed on Ubuntu 20.04 LTS using a Xeon workstation (512 GB RAM).  



---



