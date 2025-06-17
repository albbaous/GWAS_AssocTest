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

# Identify suggestive SNPs (p < 1e-5)
suggestive_snps_df <- subset(gwas, !is.na(P) & P < 1e-5)
suggestive_snps_sorted <- suggestive_snps_df[order(suggestive_snps_df$P), ]
suggestive_snps <- suggestive_snps_sorted$ID

# Print top suggestive SNPs
cat("Top Suggestive SNPs (P < 1e-5):\n")
print(suggestive_snps_sorted[, c("ID", "P")])

# Create Manhattan plot highlighting suggestive SNPs
pdf("manhattan_plot_suggestive_highlighted.pdf", width = 12, height = 6)

manhattan(
  gwas,
  chr = "CHR",
  bp = "BP",
  snp = "ID",
  p = "P",
  genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5),
  highlight = suggestive_snps,
  main = "MetaboHealth Score GWAS: Suggestive SNPs Highlighted"
)

dev.off()
