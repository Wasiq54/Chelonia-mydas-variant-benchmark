<img width="425" height="205" alt="image" src="https://github.com/user-attachments/assets/d036b7e4-e1c2-427f-928d-726e8c7ba76c" /># Methods (Tools, Versions, Commands)


##Software and Tool Versions
------------------------------------------------------------
All analyses were performed using the latest stable releases (as of 2025) in isolated Conda or Docker environments to ensure reproducibility.  
The complete list of tools and their versions is as follows:

- DeepVariant v1.9.0  
- GATK v4.6.2.0 (includes Picard v3.4.0 and HTSJDK v4.2.0)  
- FreeBayes v1.3.10  
- VarScan v2.4.6  
- BCFtools v1.22 (with HTSlib v1.22.1)  
- Samtools v1.22.1  
- BWA-MEM v0.7.17  
- Bowtie2 v2.5.2  
- Minimap2 v2.28  
- RTG Tools v3.12.1  
- FastQC v0.11.9  
- MultiQC v1.15  
- fastp v0.26.0  
- seqkit v2.8.2

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

Exact shell commands are provided in

scripts/commands_alignment.txt.

Mapping quality was evaluated using Samtools flagstat and alignment summary metrics.







Variant Calling
------------------------------------------------------------
Five variant calling tools were used to detect SNPs and indels from BWA-MEM–aligned reads:
DeepVariant v1.9.0, GATK v4.6.2.0, FreeBayes v1.3.10, VarScan v2.4.6, and BCFtools v1.22.
Each tool was executed on both normal and abnormal Chelonia mydas datasets using the same reference genome (rCheMyd1.pri.v2).
Command-line parameters for reproducibility are provided in 

scripts/commands_variant_calling.txt.

Results 

==============================================================================================================
Orignal varaints after varaint calling tools 

Tool          Sample      Total Variants
-----------------------------------------
BCFtools      Abnormal    14,380,356
BCFtools      Normal      14,327,393
DeepVariant   Abnormal    17,536,252
DeepVariant   Normal      17,425,142
FreeBayes     Abnormal    13,935,733
FreeBayes     Normal      13,910,723
GATK          Abnormal    14,671,962
GATK          Normal      14,632,668
VarScan       Abnormal    13,842,258
VarScan       Normal      13,762,684




**Variant Calling then filtering**  


Variants were retained if they:
(i) passed the caller’s internal filters (FILTER=PASS or “.”),
(ii) had alternate allele depth (AD) greater than 3, and
(iii) had variant allele frequency (VAF) above 2%.
Variants not meeting these criteria were removed.

Variants were filtered to retain only those with FILTER=PASS/., AD > 3, and VAF > 0.02.


(FILTER="PASS" or ".")
AND
ALT depth > 3
AND
VAF > 0.02)

Command-line parameters for reproducibility are provided in 
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/filter_all_tools.sh

-----------------------------------------------------------------------------        **   After Filter results **---------------------------------

Tool          Sample      Total Variants
-----------------------------------------
BCFtools      Abnormal    14,285,440
BCFtools      Normal      14,220,913
DeepVariant   Abnormal    13,767,503
DeepVariant   Normal      13,705,037
FreeBayes     Abnormal    13,332,001
FreeBayes     Normal      13,239,802
GATK          Abnormal    14,522,966
GATK          Normal      14,467,559
VarScan       Abnormal    13,833,482
VarScan       Normal      13,753,976
-----------------------------------------------------------------------------        ** SNPS files **---------------------------------
Command-line parameters for reproducibility are provided in 
(https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/make_snps.sh)


File                     SNP_File_Count
----------------------------------------
BCF_abnormal             12726145
BCF_normal               12678209
Deepvariant_abnormal     12103451
Deepvariant_normal       12058893
Freebayes_abnormal       11300566
Freebayes_normal         11247911
Gatk_abnormal            12663916
Gatk_normal              12629773
Varscan_abnormal         12288474
Varscan_normal           12228987



-----------------------------------------------------------------------------        ** Indels files **---------------------------------
Command-line parameters for reproducibility are provided in 
(https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/make_indels.sh)

File                     INDEL_File_Count
------------------------------------------
BCF_abnormal             1559295
BCF_normal               1542704
Deepvariant_abnormal     1684651
Deepvariant_normal       1666433
Freebayes_abnormal       1276508
Freebayes_normal         1254311
Gatk_abnormal            1880614
Gatk_normal              1858658
Varscan_abnormal         1545008
Varscan_normal           1524989

==================================================inter Tools concordance -----==============================
step1: make bed fils for snps and indels

Command
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/make_snp_and_indels_bed.sh


step2: SNP-only overlaps

i. Pairwise overlap (BED + bedtools intersect)    
---------------------------------------------------- normal SNPs  overlap------------------------------------------------
BCF_normal_SNPs.bed
Deepvariant_normal_SNPs.bed
Freebayes_normal_SNPs.bed
Gatk_normal_SNPs.bed
Varscan_normal_SNPs.bed


Command
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/normal_snp_overlap.sh



ToolA                         ToolB                           Overlap
----------------------------------------------------------------------------
BCF_normal_SNPs              Deepvariant_normal_SNPs          11936185
BCF_normal_SNPs              Freebayes_normal_SNPs            10988076
BCF_normal_SNPs              Gatk_normal_SNPs                 12246438
BCF_normal_SNPs              Varscan_normal_SNPs              12096398
Deepvariant_normal_SNPs      Freebayes_normal_SNPs            10566527
Deepvariant_normal_SNPs      Gatk_normal_SNPs                 11805679
Deepvariant_normal_SNPs      Varscan_normal_SNPs              11687696
Freebayes_normal_SNPs        Gatk_normal_SNPs                 10867399
Freebayes_normal_SNPs        Varscan_normal_SNPs              10648860
Gatk_normal_SNPs             Varscan_normal_SNPs              11916478


---------------------------------------------------- abnormal SNPs overlap------------------------------------------------

BCF_abnormal_SNPs.bed
Deepvariant_abnormal_SNPs.bed
Freebayes_abnormal_SNPs.bed
Gatk_abnormal_SNPs.bed
Varscan_abnormal_SNPs

Command
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/abnormal_snp_overlap.sh




ToolA                         ToolB                           Overlap
----------------------------------------------------------------------------
BCF_abnormal_SNPs            Deepvariant_abnormal_SNPs        11982196
BCF_abnormal_SNPs            Freebayes_abnormal_SNPs          11026629
BCF_abnormal_SNPs            Gatk_abnormal_SNPs               12282731
BCF_abnormal_SNPs            Varscan_abnormal_SNPs            12153920
Deepvariant_abnormal_SNPs    Freebayes_abnormal_SNPs          10607920
Deepvariant_abnormal_SNPs    Gatk_abnormal_SNPs               11849848
Deepvariant_abnormal_SNPs    Varscan_abnormal_SNPs            11744062
Freebayes_abnormal_SNPs      Gatk_abnormal_SNPs               10894738
Freebayes_abnormal_SNPs      Varscan_abnormal_SNPs            10688109
Gatk_abnormal_SNPs           Varscan_abnormal_SNPs            11959187


---------------------------------------------------- normal Indel overlap------------------------------------------------

BCF_normal_INDELs.bed
Deepvariant_normal_INDELs.bed
Freebayes_normal_INDELs.bed
Gatk_normal_INDELs.bed
Varscan_normal_INDELs.bed



Command
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/normal_indel_overlap.sh


ToolA                         ToolB                           Overlap
---------------------------------------------------------------------------
BCF_normal_INDELs            Deepvariant_normal_INDELs        1496337
BCF_normal_INDELs            Freebayes_normal_INDELs          1247269
BCF_normal_INDELs            Gatk_normal_INDELs               1500247
BCF_normal_INDELs            Varscan_normal_INDELs            1448813
Deepvariant_normal_INDELs    Freebayes_normal_INDELs          1254057
Deepvariant_normal_INDELs    Gatk_normal_INDELs               1608289
Deepvariant_normal_INDELs    Varscan_normal_INDELs            1484214
Freebayes_normal_INDELs      Gatk_normal_INDELs               1220823
Freebayes_normal_INDELs      Varscan_normal_INDELs            1189766
Gatk_normal_INDELs           Varscan_normal_INDELs            1512186




---------------------------------------------------- abnormal indel overlap------------------------------------------------
BCF_abnormal_INDELs.bed
Deepvariant_abnormal_INDELs.bed
Freebayes_abnormal_INDELs.bed
Gatk_abnormal_INDELs.bed
Varscan_abnormal_INDELs.bed



Command
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/abnormal_indel_overlap.sh


ToolA                          ToolB                            Overlap
----------------------------------------------------------------------------
BCF_abnormal_INDELs           Deepvariant_abnormal_INDELs       1513275
BCF_abnormal_INDELs           Freebayes_abnormal_INDELs         1269730
BCF_abnormal_INDELs           Gatk_abnormal_INDELs              1518467
BCF_abnormal_INDELs           Varscan_abnormal_INDELs           1467759
Deepvariant_abnormal_INDELs   Freebayes_abnormal_INDELs         1275264
Deepvariant_abnormal_INDELs   Gatk_abnormal_INDELs              1629128
Deepvariant_abnormal_INDELs   Varscan_abnormal_INDELs           1502048
Freebayes_abnormal_INDELs     Gatk_abnormal_INDELs              1239822
Freebayes_abnormal_INDELs     Varscan_abnormal_INDELs           1207299
Gatk_abnormal_INDELs          Varscan_abnormal_INDELs           1531127


----------------------------------------------------------------------------------------------------------Jaccard similarity Normal SNPs

Jaccard similarity analysis of SNP calls from the normal dataset showed a consistently high level of agreement among the five variant callers. DeepVariant, GATK, BCFtools and VarScan exhibited the strongest similarity scores, while FreeBayes showed comparatively lower but still substantial overlap. Overall, the normal‐sample SNP concordance indicates strong cross-caller consistency in high-confidence SNP regions.

command 

https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/normal_snp_jaccard.sh

----------------------------------------------------------------------------------------------------------Jaccard similarity abnormal SNPs
In the abnormal dataset, SNP-level Jaccard similarities followed a pattern similar to the normal sample. DeepVariant, GATK, BCFtools and VarScan again demonstrated strong mutual overlap, whereas FreeBayes maintained moderate similarity values. The abnormal sample displayed slightly reduced inter-caller similarity compared to the normal dataset, likely reflecting biological differences or coverage variation, but overall SNP concordance remained high

command

https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/abnormal_snp_jaccard.sh

----------------------------------------------------------------------------------------------------------Jaccard similarity Normal Indel
Jaccard similarity for INDELs in the normal dataset was lower than for SNPs, consistent with known challenges in INDEL calling. DeepVariant and GATK shared the highest INDEL concordance, followed by BCFtools and VarScan, whereas FreeBayes showed the lowest similarity with the other callers. These patterns align with previously reported INDEL-calling inconsistencies across tools.

command 
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/normal_indel_jaccard.sh


----------------------------------------------------------------------------------------------------------Jaccard similarity abnormal Indel

The abnormal dataset exhibited similar INDEL concordance patterns to the normal dataset, with DeepVariant and GATK forming the most consistent pair. Overall Jaccard values for INDELs were modest, reflecting typical inter-caller variability for small insertions and deletions. FreeBayes again showed weaker similarity to the remaining tools, while BCFtools and VarScan showed intermediate concordance.

command
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/scripts/abnormal_indel_jaccard.sh


---------------------------------------------------- For Venn SNP diagram------------------------------------------------
Venn diagrams were generated to compare the overlap of SNP calls among the five variant calling tools (DeepVariant, GATK, FreeBayes, VarScan, and BCFtools). SNPs were extracted from each filtered VCF using bcftools view, and unique variant identifiers (CHROM:POS:REF:ALT) were obtained using bcftools query. These variant sets were then imported into R, and multi-set Venn diagrams were constructed using the ggVennDiagram package to visualize shared and unique variants across tools
command 
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/venn_analysis/venn_plot_abnormal.R
https://github.com/Wasiq54/Chelonia-mydas-variant-benchmark/blob/main/venn_analysis/venn_plot_abnormal.R








### Venn Diagram Visualization of Variant Callers

To assess concordance among the five variant callers (DeepVariant, GATK, FreeBayes, VarScan, and BCFtools), SNP positions were extracted from normalized and sorted VCF files using **bcftools** and compared visually through 5-way Venn diagrams.  

The pipeline involved:
1. Normalization and indexing of VCFs using `bcftools norm`, `bcftools sort`, and `bcftools index`.
2. Extraction of SNP coordinates (`CHROM:POS`) from each caller using `bcftools query`.
3. Visualization of overlapping SNPs using the **ggVennDiagram** (R v1.3.2) package.

Scripts for complete reproducibility are available in the `venn_analysis/` directory:
- `venn_pipeline.sh` — generates normalized SNP coordinate files.  
- `venn_plot_normal.R` — visualizes overlap for the normal dataset.  
- `venn_plot_abnormal.R` — visualizes overlap for the abnormal dataset.  

The resulting figures (`FiveTool_Venn_Normal.png` and `FiveTool_Venn_Abnormal.png`) illustrate shared and unique variant calls among the tools for each dataset.


### Discordance Analysis (Missed and Unique Variants)

To evaluate discordance among variant callers, a Python-based script was developed to identify **missed** and **unique** SNPs across the five tools (DeepVariant, GATK, FreeBayes, VarScan, and BCFtools).  
The analysis used the SNP coordinate text files generated during the Venn diagram pipeline (`*_SNPs.txt`) for both normal and abnormal datasets.

**Workflow Summary**
1. Each SNP coordinate (`CHROM:POS`) was loaded from the tool-specific SNP files.
2. A combined union of all SNPs was generated.
3. For each caller:
   - **Missed SNPs** were defined as variants detected by all other tools but absent in the current one.  
   - **Unique SNPs** were defined as variants detected exclusively by that tool.
4. Counts of missed and unique SNPs were summarized per caller.

The script (`discordance_analysis.py`) is located in the `venn_analysis/` directory and can be executed as:
```bash
python venn_analysis/discordance_analysis.py
















#!/bin/bash

# ===============================
# SETTINGS
# ===============================
THREADS=60
WORKDIR="/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps"

cd "$WORKDIR"

echo "Running bcftools isec in: $WORKDIR"
echo "Threads: $THREADS"
echo "--------------------------------------"

# ===============================
# FILE LISTS
# ===============================

NORMAL_FILES=(
  BCF_normal_SNPs.vcf.gz
  Deepvariant_normal_SNPs.vcf.gz
  Freebayes_normal_SNPs.vcf.gz
  Gatk_normal_SNPs.vcf.gz
  Varscan_normal_SNPs.vcf.gz
)

ABNORMAL_FILES=(
  BCF_abnormal_SNPs.vcf.gz
  Deepvariant_abnormal_SNPs.vcf.gz
  Freebayes_abnormal_SNPs.vcf.gz
  Gatk_abnormal_SNPs.vcf.gz
  Varscan_abnormal_SNPs.vcf.gz
)

# ===============================
# FUNCTION TO RUN ISEC
# ===============================
run_isec() {
    LABEL=$1        # normal / abnormal
    N_OPTION=$2     # =3, +3, =4, +4, =5
    OUT_PREFIX=$3   # output folder name
    shift 3         

    FILES=("$@")

    echo "Running $LABEL intersection $N_OPTION → $OUT_PREFIX"

    bcftools isec \
        --threads $THREADS \
        -n "$N_OPTION" \
        -p "$OUT_PREFIX" \
        "${FILES[@]}"

    echo "Done: $OUT_PREFIX"
    echo "--------------------------------------"
}

# ===============================
# NORMAL DATASET
# ===============================

run_isec "NORMAL" "=3"  "isec_normal_n3"      "${NORMAL_FILES[@]}"
run_isec "NORMAL" "+3"  "isec_normal_n3plus"  "${NORMAL_FILES[@]}"
run_isec "NORMAL" "=4"  "isec_normal_n4"      "${NORMAL_FILES[@]}"
run_isec "NORMAL" "+4"  "isec_normal_n4plus"  "${NORMAL_FILES[@]}"
run_isec "NORMAL" "=5"  "isec_normal_n5"      "${NORMAL_FILES[@]}"

# ===============================
# ABNORMAL DATASET
# ===============================

run_isec "ABNORMAL" "=3"  "isec_abnormal_n3"      "${ABNORMAL_FILES[@]}"
run_isec "ABNORMAL" "+3"  "isec_abnormal_n3plus"  "${ABNORMAL_FILES[@]}"
run_isec "ABNORMAL" "=4"  "isec_abnormal_n4"      "${ABNORMAL_FILES[@]}"
run_isec "ABNORMAL" "+4"  "isec_abnormal_n4plus"  "${ABNORMAL_FILES[@]}"
run_isec "ABNORMAL" "=5"  "isec_abnormal_n5"      "${ABNORMAL_FILES[@]}"

echo "✨ ALL INTERSECTIONS COMPLETED SUCCESSFULLY ✨"



















===================================Truthset=====================================================

#!/usr/bin/env bash
set -euo pipefail

# Number of threads for bcftools
THREADS=40   # change to 20, 60, etc. if you want

# Base directory where your SNP VCFs live
BASE="/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps"

# Normal truthset directory (where isec_normal_* are)
NORMDIR="$BASE/truthset/Normal"

# Normal BCF VCF (used to extract full variant records)
BCF_NORMAL="$BASE/BCF_normal_SNPs.vcf.gz"

# Output directory for final truth-set VCFs
OUTDIR="$NORMDIR/truthsets"
mkdir -p "$OUTDIR"

echo "Using base VCF: $BCF_NORMAL"
echo "Writing truth sets to: $OUTDIR"
echo "Using $THREADS threads"
echo

############################
# 1) Exact n=3 consensus
############################
if [ -f "$NORMDIR/isec_normal_n3/sites.txt" ]; then
  echo "Creating truthset_normal_n3.vcf.gz (n=3 exact)..."
  bcftools view --threads $THREADS \
    -R "$NORMDIR/isec_normal_n3/sites.txt" \
    "$BCF_NORMAL" -Oz \
    -o "$OUTDIR/truthset_normal_n3.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_normal_n3.vcf.gz"
fi

############################
# 2) n>=3 (3-plus) consensus
############################
if [ -f "$NORMDIR/isec_normal_n3plus/sites.txt" ]; then
  echo "Creating truthset_normal_3plus.vcf.gz (n>=3)..."
  bcftools view --threads $THREADS \
    -R "$NORMDIR/isec_normal_n3plus/sites.txt" \
    "$BCF_NORMAL" -Oz \
    -o "$OUTDIR/truthset_normal_3plus.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_normal_3plus.vcf.gz"
fi

############################
# 3) Exact n=4 consensus
############################
if [ -f "$NORMDIR/isec_normal_n4/sites.txt" ]; then
  echo "Creating truthset_normal_n4.vcf.gz (n=4 exact)..."
  bcftools view --threads $THREADS \
    -R "$NORMDIR/isec_normal_n4/sites.txt" \
    "$BCF_NORMAL" -Oz \
    -o "$OUTDIR/truthset_normal_n4.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_normal_n4.vcf.gz"
fi

############################
# 4) n>=4 (4-plus) consensus
############################
if [ -f "$NORMDIR/isec_normal_n4plus/sites.txt" ]; then
  echo "Creating truthset_normal_4plus.vcf.gz (n>=4)..."
  bcftools view --threads $THREADS \
    -R "$NORMDIR/isec_normal_n4plus/sites.txt" \
    "$BCF_NORMAL" -Oz \
    -o "$OUTDIR/truthset_normal_4plus.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_normal_4plus.vcf.gz"
fi

############################
# 5) Exact n=5 consensus (all tools)
############################
if [ -f "$NORMDIR/isec_normal_n5/sites.txt" ]; then
  echo "Creating truthset_normal_5tools.vcf.gz (n=5, all callers)..."
  bcftools view --threads $THREADS \
    -R "$NORMDIR/isec_normal_n5/sites.txt" \
    "$BCF_NORMAL" -Oz \
    -o "$OUTDIR/truthset_normal_5tools.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_normal_5tools.vcf.gz"
fi

echo
echo "Done. Final truth sets are in: $OUTDIR"





abnoraml




#!/usr/bin/env bash
set -euo pipefail

# Number of threads for bcftools
THREADS=40   # change if you want

# Base directory where your SNP VCFs live
BASE="/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps"

# Abnormal truthset directory (where isec_abnormal_* are)
ABDIR="$BASE/truthset/abnormal"

# Abnormal BCF VCF (used to extract full variant records)
BCF_ABNORMAL="$BASE/BCF_abnormal_SNPs.vcf.gz"

# Output directory for final truth-set VCFs
OUTDIR="$ABDIR/truthsets"
mkdir -p "$OUTDIR"

echo "Using base VCF: $BCF_ABNORMAL"
echo "Writing abnormal truth sets to: $OUTDIR"
echo "Using $THREADS threads"
echo

############################
# 1) Exact n=3 consensus
############################
if [ -f "$ABDIR/isec_abnormal_n3/sites.txt" ]; then
  echo "Creating truthset_abnormal_n3.vcf.gz (n=3 exact)..."
  bcftools view --threads $THREADS \
    -R "$ABDIR/isec_abnormal_n3/sites.txt" \
    "$BCF_ABNORMAL" -Oz \
    -o "$OUTDIR/truthset_abnormal_n3.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_abnormal_n3.vcf.gz"
fi

############################
# 2) n>=3 (3-plus) consensus
############################
if [ -f "$ABDIR/isec_abnormal_n3plus/sites.txt" ]; then
  echo "Creating truthset_abnormal_3plus.vcf.gz (n>=3)..."
  bcftools view --threads $THREADS \
    -R "$ABDIR/isec_abnormal_n3plus/sites.txt" \
    "$BCF_ABNORMAL" -Oz \
    -o "$OUTDIR/truthset_abnormal_3plus.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_abnormal_3plus.vcf.gz"
fi

############################
# 3) Exact n=4 consensus
############################
if [ -f "$ABDIR/isec_abnormal_n4/sites.txt" ]; then
  echo "Creating truthset_abnormal_n4.vcf.gz (n=4 exact)..."
  bcftools view --threads $THREADS \
    -R "$ABDIR/isec_abnormal_n4/sites.txt" \
    "$BCF_ABNORMAL" -Oz \
    -o "$OUTDIR/truthset_abnormal_n4.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_abnormal_n4.vcf.gz"
fi

############################
# 4) n>=4 (4-plus) consensus
############################
if [ -f "$ABDIR/isec_abnormal_n4plus/sites.txt" ]; then
  echo "Creating truthset_abnormal_4plus.vcf.gz (n>=4)..."
  bcftools view --threads $THREADS \
    -R "$ABDIR/isec_abnormal_n4plus/sites.txt" \
    "$BCF_ABNORMAL" -Oz \
    -o "$OUTDIR/truthset_abnormal_4plus.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_abnormal_4plus.vcf.gz"
fi

############################
# 5) Exact n=5 consensus (all tools)
############################
if [ -f "$ABDIR/isec_abnormal_n5/sites.txt" ]; then
  echo "Creating truthset_abnormal_5tools.vcf.gz (n=5, all callers)..."
  bcftools view --threads $THREADS \
    -R "$ABDIR/isec_abnormal_n5/sites.txt" \
    "$BCF_ABNORMAL" -Oz \
    -o "$OUTDIR/truthset_abnormal_5tools.vcf.gz"

  bcftools index --threads $THREADS "$OUTDIR/truthset_abnormal_5tools.vcf.gz"
fi

echo
echo "Done. Abnormal truth sets are in: $OUTDIR"





















## Consensus Truth Set & Benchmarking

High-confidence variant truth sets were generated by intersecting five callers (DeepVariant, GATK, FreeBayes, VarScan, BCFtools) for normal (CH-NORMS1) and abnormal (Ab-NormS2) datasets.  

**Steps:**
1. Sort, index, and normalize VCFs using bcftools.
2. Create intersections (3+, 4+, 5+) with bcftools isec.
3. Filter caller VCFs to consensus loci using bcftools view -R.
4. Generate RTG reference SDF:
   rtg format -o reference.sdf <reference.fasta>
5. Benchmark each caller using RTG vcfeval:
   rtg vcfeval -b <truth.vcf.gz> -c <caller.vcf.gz> -t reference.sdf -o <out_dir> --vcf-score-field QUAL

All full commands are available in:

scripts/commands_consensus_truthset.txt
