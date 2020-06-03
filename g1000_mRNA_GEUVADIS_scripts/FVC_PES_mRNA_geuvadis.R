#####################

## Testing the correlation between PES and gene expression using the Geuvadis LCL RNAseq dataset

## FVC

## William Reay - April 2020 - github: https://github.com/Williamreay

#####################

set.seed(1235)

setwd("~/Desktop/Pneumonia_cytokine_lung_function/mRNA_PES_correlation/")

library(readxl)

## Load PES, PGS, and phenotype data

g1000_mRNA_FVC_PES_PGS <- read_excel("Lung_function_PES/FVC/g1000_mRNA_FVC_PES_PGS.xlsx")
g1000_mRNA_FVC_PES_PGS$SEX <- as.factor(g1000_mRNA_FVC_PES_PGS$SEX)

## Scale PES and PGS

g1000_mRNA_FVC_PES_PGS$Genome_wide_PGS_all_SNPs_threshold <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$Genome_wide_PGS_all_SNPs_threshold))
g1000_mRNA_FVC_PES_PGS$Genome_wide_PGS_0_05_threshold <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$Genome_wide_PGS_0_05_threshold))
g1000_mRNA_FVC_PES_PGS$Genome_wide_PGS_0_005_threshold <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$Genome_wide_PGS_0_005_threshold))
g1000_mRNA_FVC_PES_PGS$Pathways_in_cancer_PES <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$Pathways_in_cancer_PES))
g1000_mRNA_FVC_PES_PGS$Circadian_clock_PES <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$Circadian_clock_PES))
g1000_mRNA_FVC_PES_PGS$Class_b2_secretin_PES <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$Class_b2_secretin_PES))
g1000_mRNA_FVC_PES_PGS$NABA_matrisome_PES <- as.numeric(scale(g1000_mRNA_FVC_PES_PGS$NABA_matrisome_PES))


## Importing mRNA data ##

## Read in PEER normalised RPKM values for each individual in geuvadis

All_geuvadis_mRNA <- read.table("~/Desktop/Cross_disorder_PES_2019/GEUVADIS_gene_expression/mRNA/STABLE_ID_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt", header = TRUE)

## Merge data frames such that only European ancestry individuals with phase3 genotypes are retained

Merged_genotype_mRNA_geuvadis <- merge(All_geuvadis_mRNA, g1000_mRNA_FVC_PES_PGS, by = "ID")

## Read in the gene IDs available in the mRNA dataset

Available_ID <- read.table("~/Desktop/Cross_disorder_PES_2019/GEUVADIS_gene_expression/mRNA/Stable_IDs_GEUVADIS.txt", header = TRUE, stringsAsFactors = FALSE)

########################
## Pathways in cancer PES
########################

## Read in list of genes used to make the pathways in cancer score

Pathways_in_cancer_genes <- read.table("Preparing_geneset_specific_coordinates/Stable_IDs/Pathways_in_cancer_stable_ID.txt", header = TRUE)

## Merge pathways in cancer genes with those available in the geuvadis dataset

Merged_pathways_in_cancer <- merge(Pathways_in_cancer_genes, Available_ID, by = "Stable_ID")

List_pathways_in_cancer <- as.list(Merged_pathways_in_cancer)[[1]]

## Define linear model to test the effect of FVC Pathways in cancer PES (P < 1) on these genes

PES_pathways_in_cancer_model <- function(v){
  PES_pathways_in_cancer_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_all_SNPs_threshold + Pathways_in_cancer_PES'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PES_pathways_in_cancer_regresssion_result)) 
}

## Iterate function over list of genes

PES_pathways_in_cancer <- sapply(List_pathways_in_cancer, PES_pathways_in_cancer_model)

## Extract t-statistic and p value for each gene

PES_Extract_pathways_in_cancer <- apply(PES_pathways_in_cancer, 2, function(x) return(as.data.frame(x$coefficients)[7,3:4]))

PES_Results_pathways_in_cancer <- data.frame()
for (i in 1:length(PES_Extract_pathways_in_cancer)) {
  PES_Results_pathways_in_cancer <- rbind(PES_Results_pathways_in_cancer, PES_Extract_pathways_in_cancer[[i]])
}
rownames(PES_Results_pathways_in_cancer) <- List_pathways_in_cancer

## Apply FDR correction and export results

PES_pathways_in_cancer_FDR <- p.adjust(PES_Results_pathways_in_cancer$'Pr(>|t|)', method = "fdr")
PES_pathways_in_cancer_FINAL <- cbind(PES_Results_pathways_in_cancer, PES_pathways_in_cancer_FDR)
write.csv(PES_pathways_in_cancer_FINAL, file = "Pathways_in_cancer_PES/FVC_PES_pathways_in_cancer_all_SNPs_threshold.csv")

########################
## Circadian clock PES
########################

## Read in list of genes used to make the pathways in cancer score

Circadian_clock_genes <- read.table("Preparing_geneset_specific_coordinates/Stable_IDs/Circadian_clock_stable_IDs.txt", header = TRUE)

## Merge pathways in cancer genes with those available in the geuvadis dataset

Merged_Circadian_clock <- merge(Circadian_clock_genes, Available_ID, by = "Stable_ID")

List_Circadian_clock <- as.list(Merged_Circadian_clock)[[1]]

## Define linear model to test the effect of FVC Circadian clock PES (P < 0.05) on these genes

PES_Circadian_clock_model <- function(v){
  PES_Circadian_clock_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_0_05_threshold + Circadian_clock_PES'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PES_Circadian_clock_regresssion_result)) 
}

## Iterate function over list of genes

PES_Circadian_clock <- sapply(List_Circadian_clock, PES_Circadian_clock_model)

## Extract t-statistic and p value for each gene

PES_Extract_Circadian_clock <- apply(PES_Circadian_clock, 2, function(x) return(as.data.frame(x$coefficients)[7,3:4]))

PES_Results_Circadian_clock <- data.frame()
for (i in 1:length(PES_Extract_Circadian_clock)) {
  PES_Results_Circadian_clock <- rbind(PES_Results_Circadian_clock, PES_Extract_Circadian_clock[[i]])
}
rownames(PES_Results_Circadian_clock) <- List_Circadian_clock

## Apply FDR correction and export results

PES_Circadian_clock_FDR <- p.adjust(PES_Results_Circadian_clock$'Pr(>|t|)', method = "fdr")
PES_Circadian_clock_FINAL <- cbind(PES_Results_Circadian_clock, PES_Circadian_clock_FDR)
write.csv(PES_Circadian_clock_FINAL, file = "Circadian_clock_PES/FVC_PES_Circadian_clock_0_05_SNPs_threshold.csv")

########################
## Class b2 secretin PES
########################

## Read in list of genes used to make the pathways in cancer score

Class_b2_secretin_genes <- read.table("Preparing_geneset_specific_coordinates/Stable_IDs/Secretin_stable_IDs.txt", header = TRUE)

## Merge pathways in cancer genes with those available in the geuvadis dataset

Merged_Class_b2_secretin <- merge(Class_b2_secretin_genes, Available_ID, by = "Stable_ID")

List_Class_b2_secretin <- as.list(Merged_Class_b2_secretin)[[1]]

## Define linear model to test the effect of FVC Class b2 secretin PES (P < 0.005) on these genes

PES_Class_b2_secretin_model <- function(v){
  PES_Class_b2_secretin_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_0_005_threshold + Class_b2_secretin_PES'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PES_Class_b2_secretin_regresssion_result)) 
}

## Iterate function over list of genes

PES_Class_b2_secretin <- sapply(List_Class_b2_secretin, PES_Class_b2_secretin_model)

## Extract t-statistic and p value for each gene

PES_Extract_Class_b2_secretin <- apply(PES_Class_b2_secretin, 2, function(x) return(as.data.frame(x$coefficients)[7,3:4]))

PES_Results_Class_b2_secretin <- data.frame()
for (i in 1:length(PES_Extract_Class_b2_secretin)) {
  PES_Results_Class_b2_secretin <- rbind(PES_Results_Class_b2_secretin, PES_Extract_Class_b2_secretin[[i]])
}
rownames(PES_Results_Class_b2_secretin) <- List_Class_b2_secretin

## Apply FDR correction and export results

PES_Class_b2_secretin_FDR <- p.adjust(PES_Results_Class_b2_secretin$'Pr(>|t|)', method = "fdr")
PES_Class_b2_secretin_FINAL <- cbind(PES_Results_Class_b2_secretin, PES_Class_b2_secretin_FDR)
write.csv(PES_Class_b2_secretin_FINAL, file = "Class_b2_secretin_PES/FVC_PES_Class_b2_secretin_0_005_SNPs_threshold.csv")

########################
## NABA matrisome PES
########################

## Read in list of genes used to make the pathways in cancer score

NABA_matrisome_genes <- read.table("Preparing_geneset_specific_coordinates/Stable_IDs/NABA_matrisome_stable_IDs.txt", header = TRUE)

## Merge pathways in cancer genes with those available in the geuvadis dataset

Merged_NABA_matrisome <- merge(NABA_matrisome_genes, Available_ID, by = "Stable_ID")

List_NABA_matrisome <- as.list(Merged_NABA_matrisome)[[1]]

## Define linear model to test the effect of FVC ECM PES (P < 1) on these genes

PES_NABA_matrisome_model <- function(v){
  PES_NABA_matrisome_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_all_SNPs_threshold + NABA_matrisome_PES'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PES_NABA_matrisome_regresssion_result)) 
}

## Iterate function over list of genes

PES_NABA_matrisome <- sapply(List_NABA_matrisome, PES_NABA_matrisome_model)

## Extract t-statistic and p value for each gene

PES_Extract_NABA_matrisome <- apply(PES_NABA_matrisome, 2, function(x) return(as.data.frame(x$coefficients)[7,3:4]))

PES_Results_NABA_matrisome <- data.frame()
for (i in 1:length(PES_Extract_NABA_matrisome)) {
  PES_Results_NABA_matrisome <- rbind(PES_Results_NABA_matrisome, PES_Extract_NABA_matrisome[[i]])
}
rownames(PES_Results_NABA_matrisome) <- List_NABA_matrisome

## Apply FDR correction and export results

PES_NABA_matrisome_FDR <- p.adjust(PES_Results_NABA_matrisome$'Pr(>|t|)', method = "fdr")
PES_NABA_matrisome_FINAL <- cbind(PES_Results_NABA_matrisome, PES_NABA_matrisome_FDR)
write.csv(PES_NABA_matrisome_FINAL, file = "NABA_matrisome_PES/FVC_PES_NABA_matrisome_all_SNPs_threshold.csv")
