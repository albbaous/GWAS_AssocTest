# Load data
# Read the file
data <- read.table("final_gwas_results.MetaboHealth_Score.glm.linear", 
                   header = TRUE, sep = "\t", comment.char = "")

colnames(data)

# Subset ADD test results with valid p-values
add_data <- subset(data, TEST == "ADD" & !is.na(P) & P > 0)

# QQ plot function
qq <- function(pvector) {
  o <- -log10(sort(pvector))
  e <- -log10(ppoints(length(pvector)))
  plot(e, o, pch=19, cex=0.4, main="QQ Plot for GWAS", 
       xlab="Expected -log10(p)", ylab="Observed -log10(p)")
  abline(0, 1, col="red")
}

png("qq_plot.png", width=600, height=600)
qq(add_data$P)
dev.off()
