<img width="584" height="188" alt="image" src="https://github.com/user-attachments/assets/70836233-8317-4b29-9bca-d155653cf71d" />The raw paired-end Illumina NovaSeq 6000 reads for Chelonia mydas normal (CH-NORMS1) and abnormal (Ab-NormS2) samples were obtained from Macrogen Singapore using the TruSeq Nano DNA kit (151 bp, paired-end).
Quality evaluation with FastQC (v0.11.9) and MultiQC (v1.15) revealed no poor-quality reads, consistent base quality across all positions, and GC content near 44 %.
These results indicate high-quality DNA sequencing suitable for downstream mapping and variant analysis.


# Sequencing Data Summary

| Dataset | File | Encoding | Total Sequences | Total Bases | Read length | %GC | Poor-quality reads |
|---|---|---|---:|---:|---:|---:|---:|
| Normal (CH-NORMS1) | W1-A_1.fastq.gz | Illumina 1.9 | 456,923,725 | 68.9 Gbp | 151 bp | 44 | 0 |
| Normal (CH-NORMS1) | W1-A_2.fastq.gz | Illumina 1.9 | 456,923,725 | 68.9 Gbp | 151 bp | 45 | 0 |
| Abnormal (Ab-NormS2) | S1_1.fastq.gz | Illumina 1.9 | 441,497,634 | 66.6 Gbp | 151 bp | 44 | 0 |
| Abnormal (Ab-NormS2) | S1_2.fastq.gz | Illumina 1.9 | 441,497,634 | 66.6 Gbp | 151 bp | 44 | 0 |






File                                  Total        SNPs         INDELs
---------------------------------------------------------------------------
BCF_abnormal.vcf.gz                   14380356     12801610     1578746
BCF_normal.vcf.gz                     14327393     12760468     1566925
Deepvariant_abnormal.vcf.gz           17536252     15252776     2304353
Deepvariant_normal.vcf.gz             17425142     15191513     2254495
Freebayes_abnormal.vcf.gz            13935733     11776492     1322191
Freebayes_normal.vcf.gz              13910723     11792036     1307234
Gatk_abnormal.vcf.gz                 14671962     12774695     1919568
Gatk_normal.vcf.gz                   14632668     12749971     1904422
Varscan_abnormal.vcf.gz              13842258     12296877     1545381
Varscan_normal.vcf.gz                13762684     12237274     1525410

               After Filter 


Tool        Sample     File                                   Total        SNPs         INDELs
-------------------------------------------------------------------------------------------------------
BCFtools    Abnormal   BCF_abnormal.filtered.vcf.gz           14285440     12726145     1559295
BCFtools    Normal     BCF_normal.filtered.vcf.gz             14220913     12678209     1542704
DeepVariant Abnormal   Deepvariant_abnormal.filtered.vcf.gz   13767503     12103451     1684651
DeepVariant Normal     Deepvariant_normal.filtered.vcf.gz     13705037     12058893     1666433
FreeBayes   Abnormal   Freebayes_abnormal.filtered.vcf.gz     13332001     11300566     1276508
FreeBayes   Normal     Freebayes_normal.filtered.vcf.gz       13239802     11247911     1254311
GATK        Abnormal   Gatk_abnormal.filtered.vcf.gz          14522966     12663916     1880614
GATK        Normal     Gatk_normal.filtered.vcf.gz            14467559     12629773     1858658
VarScan     Abnormal   Varscan_abnormal.filtered.vcf.gz       13833482     12288474     1545008
VarScan     Normal     Varscan_normal.filtered.vcf.gz         13753976     12228987     1524989


