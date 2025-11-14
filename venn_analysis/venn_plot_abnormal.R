#!/usr/bin/env Rscript
# Five-caller Venn (ABNORMAL dataset)
# Requires: ggVennDiagram, ggplot2, svglite
# install.packages(c("ggVennDiagram", "ggplot2", "svglite"))

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
})

# Path to abnormal SNP BED directory
abnormal_dir <- "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed"
setwd(abnormal_dir)

# Helper function: read BED and convert to "chr:pos"
read_bed_ids <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)
  unique(paste(df$V1, df$V2, sep = ":"))
}

# Read ABNORMAL SNP files
venn_list <- list(
  BCF         = read_bed_ids("BCF_abnormal_SNPs.bed"),
  DeepVariant = read_bed_ids("Deepvariant_abnormal_SNPs.bed"),
  FreeBayes   = read_bed_ids("Freebayes_abnormal_SNPs.bed"),
  GATK        = read_bed_ids("Gatk_abnormal_SNPs.bed"),
  VarScan     = read_bed_ids("Varscan_abnormal_SNPs.bed")
)

# Print counts for log
cat("Total SNPs per caller (ABNORMAL):\n")
print(sapply(venn_list, length))

# Generate clean non-overlapping Venn diagram
p <- ggVennDiagram(
  venn_list,
  label_alpha = 0,      # solid label background
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  ggtitle("Five-Tool Venn Diagram (Abnormal SNPs)")

# Save PNG (high resolution)
ggsave(
  filename = "FiveTool_Venn_Abnormal.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 600
)

# Save SVG (vector format, perfect for publication)
ggsave(
  filename = "FiveTool_Venn_Abnormal.svg",
  plot = p,
  width = 10,
  height = 8,
  dpi = 600
)

cat("Saved: FiveTool_Venn_Abnormal.png & FiveTool_Venn_Abnormal.svg in ", abnormal_dir, "\n", sep = "")

