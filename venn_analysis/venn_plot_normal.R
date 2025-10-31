#!/usr/bin/env Rscript
# Five-caller Venn (NORMAL dataset)
# Requires: ggVennDiagram, ggplot2, svglite
# install.packages(c("ggVennDiagram", "ggplot2", "svglite"))

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
  library(svglite)
})

normal_dir <- "/media/fahad/DNA_Work2/concordance/ForvanDiagram"
setwd(normal_dir)

read_snps <- function(f) if (file.exists(f)) scan(f, what = character(), quiet = TRUE) else character()

snps <- list(
  BCFTools    = read_snps("BCFTools_Normal_SNPs.txt"),
  DeepVariant = read_snps("DeepVariant_Normal_SNPs.txt"),
  FreeBayes   = read_snps("FreeBayes_Normal_SNPs.txt"),
  GATK        = read_snps("GATK_Normal_SNPs.txt"),
  VarScan     = read_snps("VarScan_Normal_SNPs.txt")
)

p <- ggVennDiagram(snps, label_alpha = 0, label = "count") +
  theme_void() + theme(legend.position = "none")

ggsave("FiveTool_Venn_Normal.png", p, width = 10, height = 8, dpi = 600)
ggsave("FiveTool_Venn_Normal.svg", p, width = 10, height = 8, dpi = 600)
cat("Saved: FiveTool_Venn_Normal.png/.svg in ", normal_dir, "\n", sep = "")
