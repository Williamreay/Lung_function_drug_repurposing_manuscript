#####################

## Correlation plot of PES for the FEV1 and FVC gene-sets considered

## William Reay - May 2020 - github: https://github.com/Williamreay

#####################
set.seed(1235)

library(readxl)
library(corrplot)

setwd("~/Desktop/Pneumonia_cytokine_lung_function/mRNA_PES_correlation/")

## FEV1

FEV1_corr_plot <- read_excel("Lung_function_PES/FEV1_corr_plot.xlsx")

## Scale PES and PGS

FEV1_corr_plot$Genome_wide_PGS_all_SNPs <- as.numeric(scale(FEV1_corr_plot$Genome_wide_PGS_all_SNPs))
FEV1_corr_plot$Genome_wide_PGS_0.05 <- as.numeric(scale(FEV1_corr_plot$Genome_wide_PGS_0.05))
FEV1_corr_plot$Dilated_cardiomyopathy <- as.numeric(scale(FEV1_corr_plot$Dilated_cardiomyopathy))
FEV1_corr_plot$Pathways_in_cancer <- as.numeric(scale(FEV1_corr_plot$Pathways_in_cancer))
FEV1_corr_plot$Extension_of_telomeres <- as.numeric(scale(FEV1_corr_plot$Extension_of_telomeres))

## Generate correlation plot

FEV1_Corr_plot_input <- FEV1_corr_plot

All_PES_Corr <- cor(FEV1_Corr_plot_input)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(All_PES_Corr, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.45, number.cex= 0.5)

## Generate correlation plot which only shows significant correlations

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(FEV1_Corr_plot_input)

corrplot(All_PES_Corr, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.6, number.cex= 0.5, p.mat = p.mat, sig.level = 0.05, insig = "blank")

## FVC

FVC_corr_plot <- read_excel("Lung_function_PES/FVC_corr_plot.xlsx")

## Scale PES and PGS

FVC_corr_plot$Genome_wide_PGS_all_SNPs <- as.numeric(scale(FVC_corr_plot$Genome_wide_PGS_all_SNPs))
FVC_corr_plot$Genome_wide_PGS_0.05 <- as.numeric(scale(FVC_corr_plot$Genome_wide_PGS_0.05))
FVC_corr_plot$Genome_wide_PGS_0.005 <- as.numeric(scale(FVC_corr_plot$Genome_wide_PGS_0.005))
FVC_corr_plot$Circadian_clock <- as.numeric(scale(FVC_corr_plot$Circadian_clock))
FVC_corr_plot$Pathways_in_cancer <- as.numeric(scale(FVC_corr_plot$Pathways_in_cancer))
FVC_corr_plot$Extracellular_matrix <- as.numeric(scale(FVC_corr_plot$Extracellular_matrix))
FVC_corr_plot$Class_b2_secretin <- as.numeric(scale(FVC_corr_plot$Class_b2_secretin))

FVC_Corr_plot_input <- FVC_corr_plot

All_PES_Corr <- cor(FVC_Corr_plot_input)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

p.mat <- cor.mtest(FVC_Corr_plot_input)

corrplot(All_PES_Corr, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.6, number.cex= 0.5, p.mat = p.mat, sig.level = 0.00102040816, insig = "blank")
