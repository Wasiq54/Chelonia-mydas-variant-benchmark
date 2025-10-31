#!/usr/bin/env Rscript
# Five-caller Venn (ABNORMAL dataset)
# Requires: ggVennDiagram, ggplot2, svglite
# install.packages(c("ggVennDiagram", "ggplot2", "svglite"))

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
  library(svglite)
})

ab_dir <- "/media/fahad/DNA_Work2/concordance/ForvanDiagram/Abnormal"
setwd(ab_dir)

read_snps <- function(f) if (file.exists(f)) scan(f, what = character(), quiet = TRUE) else character()

snps <- list(
  BCFTools    = read_snps("BCFTools_Abnorm_SNPs.txt"),
  DeepVariant = read_snps("DeepVariant_Abnorm_SNPs.txt"),
  FreeBayes   = read_snps("FreeBayes_Abnorm_SNPs.txt"),
  GATK        = read_snps("GATK_Abnorm_SNPs.txt"),
  VarScan     = read_snps("VarScan_Abnorm_SNPs.txt")
)

p <- ggVennDiagram(snps, label_alpha = 0, label = "count") +
  theme_void() + theme(legend.position = "none")

ggsave("FiveTool_Venn_Abnormal.png", p, width = 10, height = 8, dpi = 600)
ggsave("FiveTool_Venn_Abnormal.svg", p, width = 10, height = 8, dpi = 600)
cat("Saved: FiveTool_Venn_Abnormal.png/.svg in ", ab_dir, "\n", sep = "")
