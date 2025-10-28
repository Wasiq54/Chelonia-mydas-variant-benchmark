# Chelonia-mydas-variant-benchmark
Benchmarking DeepVariant and conventional variant callers for the Chelonia mydas genome — including scripts, workflows, and reproducibility materials.


This repo contains all steps & scripts for:
“Benchmarking DeepVariant and Conventional Variant Callers for Consensus-Based Analysis of the Chelonia mydas Genome.”

## Data access (to be updated)
- ENA BioProject: PRJNAxxxxxx
- EVA study (VCFs): PRJEByyyyy
- Archived code (Zenodo DOI): 10.5281/zenodo.zzzzzz

## Workflow (high level)
1) QC: FastQC + MultiQC
2) (Optional) Trim: fastp
3) Align: BWA-MEM (and Bowtie2/Minimap2 for comparison)
4) Call variants: DeepVariant, GATK, FreeBayes, VarScan, BCFtools
5) Evaluate: RTG vcfeval (+ tables & plots)

See **methods.md** for exact commands & versions.

