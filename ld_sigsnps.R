library(dplyr)

# Example: GWAS data frame should have columns: CHR, BP, P, ID
# gwas <- read.table("your_file.txt", header=TRUE)

# Define your regions as a named list: c(start, end)
regions <- list(
  region1 = c(21560000, 21650000),  # Main cluster
  region2 = c(35600000, 35700000),  # Second cluster
  region3 = c(49980000, 50000000),  # Third cluster
  region4 = c(19970000, 20000000)   # Fourth cluster
)

# Extract lead SNP (smallest P) from each region
lead_snps <- lapply(regions, function(region) {
  gwas %>%
    filter(BP >= region[1], BP <= region[2]) %>%
    arrange(P) %>%
    slice(1)
}) %>% bind_rows()

# View or output the lead SNPs
print(lead_snps)
write.csv(lead_snps, "lead_snps_table.csv", row.names=FALSE)

# Assume your data frame is called 'gwas' and has columns: ID, CHR, BP, P

# Find the SNP with the smallest p-value
lead_snp <- gwas[which.min(gwas$P), ]

print(lead_snp)

