#!/usr/bin/env Rscript
# Five-caller Venn (NORMAL dataset)
# Requires: ggVennDiagram, ggplot2, svglite
# install.packages(c("ggVennDiagram", "ggplot2", "svglite"))

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
})

setwd("/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed")

read_bed_ids <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)
  unique(paste(df$V1, df$V2, sep=":"))
}

venn_list <- list(
  BCF         = read_bed_ids("BCF_normal_SNPs.bed"),
  DeepVariant = read_bed_ids("Deepvariant_normal_SNPs.bed"),
  FreeBayes   = read_bed_ids("Freebayes_normal_SNPs.bed"),
  GATK        = read_bed_ids("Gatk_normal_SNPs.bed"),
  VarScan     = read_bed_ids("Varscan_normal_SNPs.bed")
)

cat("Total SNPs per caller (NORMAL):\n")
print(sapply(venn_list, length))

p <- ggVennDiagram(
  venn_list,
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none")

ggsave(
  "normal_SNPs_5callers_VENN.png",
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)

cat("Saved: normal_SNPs_5callers_VENN.png\n")

