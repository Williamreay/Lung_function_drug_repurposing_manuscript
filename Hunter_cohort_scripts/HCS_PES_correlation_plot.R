#############################
## Correlation plot of FEV1 and FVC PES in the Hunter Community Study Cohort

## William Reay - May 2020 - github: https://github.com/Williamreay
#############################

library(readxl)
library(corrplot)

setwd("~/Desktop/Pneumonia_cytokine_lung_function/HCS_PES_profiles/")

## FEV1

FEV1_corr_plot <- read_excel("FEV1_HCS_PES_combined.xlsx")

## Scale PES and PGS
FEV1_corr_plot$Genome_wide_all_SNPs_threshold <- as.numeric(scale(FEV1_corr_plot$Genome_wide_all_SNPs_threshold))
FEV1_corr_plot$Genome_wide_0.05_threshold <- as.numeric(scale(FEV1_corr_plot$Genome_wide_0.05_threshold))
FEV1_corr_plot$Dilated_cardiomyopathy_PES <- as.numeric(scale(FEV1_corr_plot$Dilated_cardiomyopathy_PES))
FEV1_corr_plot$Pathways_in_cancer_PES <- as.numeric(scale(FEV1_corr_plot$Pathways_in_cancer_PES))
FEV1_corr_plot$Extension_of_telomeres_PES <- as.numeric(scale(FEV1_corr_plot$Extension_of_telomeres_PES))

## Generate correlation plot

FEV1_Corr_plot_input <- FEV1_corr_plot[-c(1:6)]

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
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.6, number.cex= 0.5, p.mat = p.mat, sig.level = 0.002, insig = "blank")


## FVC

FVC_corr_plot <- read_excel("FVC_HCS_PES_combined.xlsx")

## Scale PES and PGS

FVC_corr_plot$Genome_wide_all_SNPs_threshold <- as.numeric(scale(FVC_corr_plot$Genome_wide_all_SNPs_threshold))
FVC_corr_plot$Genome_wide_0.05_threshold <- as.numeric(scale(FVC_corr_plot$Genome_wide_0.05_threshold))
FVC_corr_plot$Genome_wide_0.005_threshold <- as.numeric(scale(FVC_corr_plot$Genome_wide_0.005_threshold))
FVC_corr_plot$Circadian_clock_PES <- as.numeric(scale(FVC_corr_plot$Circadian_clock_PES))
FVC_corr_plot$Pathways_in_cancer_PES <- as.numeric(scale(FVC_corr_plot$Pathways_in_cancer_PES))
FVC_corr_plot$ECM_PES <- as.numeric(scale(FVC_corr_plot$ECM_PES))
FVC_corr_plot$Class_b2_Secretin_PES <- as.numeric(scale(FVC_corr_plot$Class_b2_Secretin_PES))

FVC_Corr_plot_input <- FVC_corr_plot[-c(1:6)]

All_PES_Corr <- cor(FVC_Corr_plot_input)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

p.mat <- cor.mtest(FVC_Corr_plot_input)

corrplot(All_PES_Corr, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.6, number.cex= 0.5, p.mat = p.mat, sig.level = 0.00102040816, insig = "blank")

##################
## Univariate regression and VIF for each PES - PGS correlation (adjusted for first 5 PCs)
##################

## FEV1

## Pathways in cancer vs FEV1 PGS (all SNPs)

Cancer_FEV1 <- lm(Pathways_in_cancer_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, data = FEV1_corr_plot)
summary(Cancer_FEV1)
car::vif(Cancer_FEV1)

## Dilated cardiomyopathy vs FEV1 PGS (all SNPS)

Dilated_cardiomyopathy_FEV1 <- lm(Dilated_cardiomyopathy_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, data = FEV1_corr_plot)
summary(Dilated_cardiomyopathy_FEV1)
car::vif(Dilated_cardiomyopathy_FEV1)

## Extension of telomeres vs FEV1 PGS (P < 0.05)
Telomeres_FEV1 <- lm(Extension_of_telomeres_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.05_threshold, data = FEV1_corr_plot)
summary(Telomeres_FEV1)
car::vif(Telomeres_FEV1)

## FVC

## Pathways in cancer vs FVC PGS (all SNPs)
Cancer_FVC <- lm(Pathways_in_cancer_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, data = FVC_corr_plot)
summary(Cancer_FVC)
car::vif(Cancer_FVC)

## Circadian clock vs FVC PGS (P < 0.05)
Circadian_FVC <- lm(Circadian_clock_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.05_threshold, data = FVC_corr_plot)
summary(Circadian_FVC)
car::vif(Circadian_FVC)

## Class b2 Secretin vs FVC PGS (P < 0.005)
Secretin_FVC <- lm(Class_b2_Secretin_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.005_threshold, data = FVC_corr_plot)
summary(Secretin_FVC)
car::vif(Secretin_FVC)

## Extracellular matrix (ECM) vs FVC PGS (all SNPs)
ECM_FVC <- lm(ECM_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, data = FVC_corr_plot)
summary(ECM_FVC)
car::vif(ECM_FVC)
