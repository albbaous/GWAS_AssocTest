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

# Calculate Benjamini-Hochberg (FDR) adjusted p-values
gwas$FDR_BH <- p.adjust(gwas$P, method = "BH")

# Set your FDR significance threshold (commonly 0.05)
fdr_threshold <- 0.05

# Identify FDR-significant SNPs
sig_fdr_df    <- subset(gwas, !is.na(FDR_BH) & FDR_BH < fdr_threshold)
sig_fdr_sorted <- sig_fdr_df[order(sig_fdr_df$FDR_BH), ]

# Count the number of FDR-significant SNPs
num_fdr_snps <- nrow(sig_fdr_sorted)
cat("Number of SNPs identified as FDR-significant (FDR <", fdr_threshold, "):", num_fdr_snps, "\n")

cat("Benjamini-Hochberg FDR-significant SNPs (FDR <", fdr_threshold, "):\n")
print(sig_fdr_sorted[, c("ID", "CHR", "BP", "P", "FDR_BH")])

# Write all FDR-significant SNPs to a tab-delimited file
write.table(
  sig_fdr_sorted[, c("ID", "CHR", "BP", "P", "FDR_BH")],
  file = "FDR_significant_SNPs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Create Manhattan plot highlighting FDR-significant SNPs
pdf("manhattan_plot_fdr_highlighted.pdf", width = 12, height = 6)

manhattan(
  gwas,
  chr           = "CHR",
  bp            = "BP",
  snp           = "ID",
  p             = "P",
  genomewideline = FALSE,  # No Bonferroni line here, but can add a custom -log10(FDR threshold) if you wish
  suggestiveline = FALSE,
  highlight     = sig_fdr_sorted$ID,
  main          = "MetaboHealth Score GWAS: FDR-significant SNPs"
)

dev.off()

