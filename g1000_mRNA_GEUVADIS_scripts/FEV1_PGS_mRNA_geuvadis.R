#####################

## Testing the correlation between genome-wide PGS and gene expression using the Geuvadis LCL RNAseq dataset

## FEV1

## William Reay - April 2020 - github: https://github.com/Williamreay

#####################

set.seed(1235)

setwd("~/Desktop/Pneumonia_cytokine_lung_function/mRNA_PES_correlation/")

library(readxl)

## Load PES, PGS, and phenotype data

g1000_mRNA_FEV1_PES_PGS <- read_excel("Lung_function_PES/FEV1/g1000_mRNA_FEV1_PES_PGS.xlsx")
g1000_mRNA_FEV1_PES_PGS$SEX <- as.factor(g1000_mRNA_FEV1_PES_PGS$SEX)

## Scale PES and PGS

g1000_mRNA_FEV1_PES_PGS$Genome_wide_PGS_all_SNPs_threshold <- as.numeric(scale(g1000_mRNA_FEV1_PES_PGS$Genome_wide_PGS_all_SNPs_threshold))
g1000_mRNA_FEV1_PES_PGS$Genome_wide_PGS_0_5_threshold <- as.numeric(scale(g1000_mRNA_FEV1_PES_PGS$Genome_wide_PGS_0_5_threshold))
g1000_mRNA_FEV1_PES_PGS$Genome_wide_PGS_0_05_threshold <- as.numeric(scale(g1000_mRNA_FEV1_PES_PGS$Genome_wide_PGS_0_05_threshold))

## Importing mRNA data ##

## Read in PEER normalised RPKM values for each individual in geuvadis

All_geuvadis_mRNA <- read.table("~/Desktop/Cross_disorder_PES_2019/GEUVADIS_gene_expression/mRNA/STABLE_ID_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt", header = TRUE)

## Merge data frames such that only European ancestry individuals with phase3 genotypes are retained

Merged_genotype_mRNA_geuvadis <- merge(All_geuvadis_mRNA, g1000_mRNA_FEV1_PES_PGS, by = "ID")

## Read in the gene IDs available in the mRNA dataset

Available_ID <- read.table("~/Desktop/Cross_disorder_PES_2019/GEUVADIS_gene_expression/mRNA/Stable_IDs_GEUVADIS.txt", header = TRUE, stringsAsFactors = FALSE)

########################
## Genome wide PGS - all SNPs
#########################

List_available_ID <- as.list(Available_ID)[[1]]

## Define linear model to test the effect of FEV1 PGS (P < 1) on these genes

PGS_all_SNPs_model <- function(v){
  PGS_all_SNPs_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_all_SNPs_threshold'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PGS_all_SNPs_regresssion_result)) 
}

## Iterate function over list of genes

PGS_all_SNPs <- sapply(List_available_ID, PGS_all_SNPs_model)

## Extract t-statistic and p value for each gene

PGS_Extract_all_SNPs <- apply(PGS_all_SNPs, 2, function(x) return(as.data.frame(x$coefficients)[6,3:4]))

PGS_Results_all_SNPs <- data.frame()
for (i in 1:length(PGS_Extract_all_SNPs)) {
  PGS_Results_all_SNPs <- rbind(PGS_Results_all_SNPs, PGS_Extract_all_SNPs[[i]])
}
rownames(PGS_Results_all_SNPs) <- List_available_ID

## Apply FDR correction and export results

PGS_all_SNPs_FDR <- p.adjust(PGS_Results_all_SNPs$'Pr(>|t|)', method = "fdr")

PGS_all_SNPs_FINAL <- cbind(PGS_Results_all_SNPs, PGS_all_SNPs_FDR)
write.csv(PGS_all_SNPs_FINAL, file = "Genome_wide_PGS_results/FEV1_PGS_all_SNPs.csv")

########################
## Genome wide PGS - P < 0.05
#########################

## Define linear model to test the effect of FEV1 PGS (P < 0.05) on these genes

PGS_0_05_threshold_model <- function(v){
  PGS_0_05_threshold_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_0_05_threshold'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PGS_0_05_threshold_regresssion_result)) 
}

## Iterate function over list of genes

PGS_0_05_threshold <- sapply(List_available_ID, PGS_0_05_threshold_model)

## Extract t-statistic and p value for each gene

PGS_Extract_0_05_threshold <- apply(PGS_0_05_threshold, 2, function(x) return(as.data.frame(x$coefficients)[6,3:4]))

PGS_Results_0_05_threshold <- data.frame()
for (i in 1:length(PGS_Extract_0_05_threshold)) {
  PGS_Results_0_05_threshold <- rbind(PGS_Results_0_05_threshold, PGS_Extract_0_05_threshold[[i]])
}
rownames(PGS_Results_0_05_threshold) <- List_available_ID

## Apply FDR correction and export results

PGS_0_05_threshold_FDR <- p.adjust(PGS_Results_0_05_threshold$'Pr(>|t|)', method = "fdr")
PGS_0_05_threshold_FINAL <- cbind(PGS_Results_0_05_threshold, PGS_0_05_threshold_FDR)
write.csv(PGS_0_05_threshold_FINAL, file = "Genome_wide_PGS_results/FEV1_PGS_0_05_threshold.csv")

remove(list = ls())
