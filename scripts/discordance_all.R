


suppressPackageStartupMessages({
  library(tools)
})

# --------- Helper: read BED and make ID like "chr:pos" ----------
read_bed_ids <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)
  unique(paste(df$V1, df$V2, sep=":"))
}

# --------- Define your datasets (edit paths if needed) ----------
datasets <- list(

  normal_SNPs = list(
    files = c(
      BCF         = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/BCF_normal_SNPs.bed",
      DeepVariant = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Deepvariant_normal_SNPs.bed",
      FreeBayes   = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Freebayes_normal_SNPs.bed",
      GATK        = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Gatk_normal_SNPs.bed",
      VarScan     = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Varscan_normal_SNPs.bed"
    ),
    out = "discordance_normal_SNPs.tsv"
  ),

  abnormal_SNPs = list(
    files = c(
      BCF         = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/BCF_abnormal_SNPs.bed",
      DeepVariant = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Deepvariant_abnormal_SNPs.bed",
      FreeBayes   = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Freebayes_abnormal_SNPs.bed",
      GATK        = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Gatk_abnormal_SNPs.bed",
      VarScan     = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/snps_bed/Varscan_abnormal_SNPs.bed"
    ),
    out = "discordance_abnormal_SNPs.tsv"
  ),

  normal_INDELs = list(
    files = c(
      BCF         = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/BCF_normal_INDELs.bed",
      DeepVariant = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Deepvariant_normal_INDELs.bed",
      FreeBayes   = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Freebayes_normal_INDELs.bed",
      GATK        = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Gatk_normal_INDELs.bed",
      VarScan     = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Varscan_normal_INDELs.bed"
    ),
    out = "discordance_normal_INDELs.tsv"
  ),

  abnormal_INDELs = list(
    files = c(
      BCF         = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/BCF_abnormal_INDELs.bed",
      DeepVariant = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Deepvariant_abnormal_INDELs.bed",
      FreeBayes   = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Freebayes_abnormal_INDELs.bed",
      GATK        = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Gatk_abnormal_INDELs.bed",
      VarScan     = "/home/work/Desktop/variants/new_latest_all_data/snvs/working/indels_bed/Varscan_abnormal_INDELs.bed"
    ),
    out = "discordance_abnormal_INDELs.tsv"
  )
)

# --------- Main loop: compute discordance for each dataset ----------
for (ds in names(datasets)) {
  info <- datasets[[ds]]

  # read each tool as a set of IDs
  venn_list <- lapply(info$files, read_bed_ids)
  names(venn_list) <- names(info$files)

  # total variants per tool
  total <- sapply(venn_list, length)

  # variants shared by ALL 5 tools (fully concordant core)
  all5_ids <- Reduce(intersect, venn_list)
  shared_all5 <- length(all5_ids)

  # unique variants: present ONLY in this tool (FP-type discordance)
  unique_counts <- sapply(names(venn_list), function(nm) {
    this   <- venn_list[[nm]]
    others <- Reduce(union, venn_list[names(venn_list) != nm])
    length(setdiff(this, others))
  })

  # missed vs ALL-5 core: variants in all5 but absent in this tool
  missed_vs_all5 <- sapply(names(venn_list), function(nm) {
    length(setdiff(all5_ids, venn_list[[nm]]))
  })

  # missed vs UNION of other tools (broader FN-type discordance)
  missed_vs_union_others <- sapply(names(venn_list), function(nm) {
    this         <- venn_list[[nm]]
    others_union <- Reduce(union, venn_list[names(venn_list) != nm])
    length(setdiff(others_union, this))
  })

  # build summary table
  out_df <- data.frame(
    Dataset                 = ds,
    Tool                    = names(venn_list),
    Total                   = as.integer(total),
    Shared_all5             = as.integer(rep(shared_all5, length(total))),
    Unique_only_this        = as.integer(unique_counts),
    Missed_vs_all5          = as.integer(missed_vs_all5),
    Missed_vs_union_others  = as.integer(missed_vs_union_others),
    stringsAsFactors = FALSE
  )

  # write to file in current working dir
  write.table(out_df,
              file = info$out,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  cat("Saved:", info$out, "\n")
}
