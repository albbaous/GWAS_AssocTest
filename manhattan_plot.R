# Set working directory
setwd("/Users/user/Desktop/Biostatistics1/UKB visualisation/BGEN/gwas")

# Load required package
library(qqman)

# Read GWAS results
gwas <- read.table(
  "final_gwas_results.MetaboHealth_Score.glm.linear",
  header = TRUE,
  sep = "\t",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Rename chromosome column
names(gwas)[names(gwas) == "#CHROM"] <- "CHR"

# Filter to additive model only
gwas <- subset(gwas, TEST == "ADD")

# Clean/convert necessary columns
gwas$CHR <- as.numeric(gwas$CHR)
gwas$BP  <- as.numeric(gwas$POS)
gwas$P   <- as.numeric(gwas$P)
gwas$ID  <- trimws(gwas$ID)

# Calculate Bonferroni threshold for 570,731 SNPs
n_snps <- 570731
bonf_threshold <- 0.05 / n_snps
cat("Bonferroni significance threshold:", format(bonf_threshold, scientific = TRUE), "\n")

# Identify Bonferroni-significant SNPs
sig_snps_df    <- subset(gwas, !is.na(P) & P < bonf_threshold)
sig_snps_sorted <- sig_snps_df[order(sig_snps_df$P), ]

# Print Bonferroni-significant SNPs
cat("Bonferroni-significant SNPs (P <", format(bonf_threshold, scientific = TRUE), "):\n")
print(sig_snps_sorted[, c("ID", "CHR", "BP", "P")])

# Create Manhattan plot highlighting only Bonferroni-significant SNPs
pdf("manhattan_plot_bonferroni_highlighted.pdf", width = 12, height = 6)

manhattan(
  gwas,
  chr           = "CHR",
  bp            = "BP",
  snp           = "ID",
  p             = "P",
  genomewideline = -log10(bonf_threshold),
  suggestiveline = FALSE,
  highlight     = sig_snps_sorted$ID,
  main          = "MetaboHealth Score GWAS: Bonferroni-significant SNPs"
)

dev.off()
