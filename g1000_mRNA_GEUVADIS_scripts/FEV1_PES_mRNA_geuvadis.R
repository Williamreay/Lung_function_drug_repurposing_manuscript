#####################

## Testing the correlation between PES and gene expression using the Geuvadis LCL RNAseq dataset

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
g1000_mRNA_FEV1_PES_PGS$Dilated_cardiomyopathy_PES <- as.numeric(scale(g1000_mRNA_FEV1_PES_PGS$Dilated_cardiomyopathy_PES))
g1000_mRNA_FEV1_PES_PGS$Pathways_in_cancer_PES <- as.numeric(scale(g1000_mRNA_FEV1_PES_PGS$Pathways_in_cancer_PES))
g1000_mRNA_FEV1_PES_PGS$Extension_of_telomeres_PES <- as.numeric(scale(g1000_mRNA_FEV1_PES_PGS$Extension_of_telomeres_PES))


## Importing mRNA data ##

## Read in PEER normalised RPKM values for each individual in geuvadis

All_geuvadis_mRNA <- read.table("~/Desktop/Cross_disorder_PES_2019/GEUVADIS_gene_expression/mRNA/STABLE_ID_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt", header = TRUE)

## Merge data frames such that only European ancestry individuals with phase3 genotypes are retained

Merged_genotype_mRNA_geuvadis <- merge(All_geuvadis_mRNA, g1000_mRNA_FEV1_PES_PGS, by = "ID")

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

## Define linear model to test the effect of FEV1 Pathways in cancer PES (P < 1) on these genes

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
write.csv(PES_pathways_in_cancer_FINAL, file = "Pathways_in_cancer_PES/FEV1_PES_pathways_in_cancer_all_SNPs_threshold.csv")

########################
## Dilated cardiomyopathy PES
########################

## Read in list of genes used to make the pathways in cancer score

Dilated_cardiomyopathy_genes <- read.table("Preparing_geneset_specific_coordinates/Stable_IDs/Dilated_cardiomyopathy_stable_IDs.txt", header = TRUE)

## Merge pathways in cancer genes with those available in the geuvadis dataset

Merged_Dilated_cardiomyopathy <- merge(Dilated_cardiomyopathy_genes, Available_ID, by = "Stable_ID")

List_Dilated_cardiomyopathy <- as.list(Merged_Dilated_cardiomyopathy)[[1]]

## Define linear model to test the effect of FEV1 Pathways in cancer PES (P < 1) on these genes

PES_Dilated_cardiomyopathy_model <- function(v){
  PES_Dilated_cardiomyopathy_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_all_SNPs_threshold + Dilated_cardiomyopathy_PES'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PES_Dilated_cardiomyopathy_regresssion_result)) 
}

## Iterate function over list of genes

PES_Dilated_cardiomyopathy <- sapply(List_Dilated_cardiomyopathy, PES_Dilated_cardiomyopathy_model)

## Extract t-statistic and p value for each gene

PES_Extract_Dilated_cardiomyopathy <- apply(PES_Dilated_cardiomyopathy, 2, function(x) return(as.data.frame(x$coefficients)[7,3:4]))

PES_Results_Dilated_cardiomyopathy <- data.frame()
for (i in 1:length(PES_Extract_Dilated_cardiomyopathy)) {
  PES_Results_Dilated_cardiomyopathy <- rbind(PES_Results_Dilated_cardiomyopathy, PES_Extract_Dilated_cardiomyopathy[[i]])
}
rownames(PES_Results_Dilated_cardiomyopathy) <- List_Dilated_cardiomyopathy

## Apply FDR correction and export results

PES_Dilated_cardiomyopathy_FDR <- p.adjust(PES_Results_Dilated_cardiomyopathy$'Pr(>|t|)', method = "fdr")
PES_Dilated_cardiomyopathy_FINAL <- cbind(PES_Results_Dilated_cardiomyopathy, PES_Dilated_cardiomyopathy_FDR)
write.csv(PES_Dilated_cardiomyopathy_FINAL, file = "Dilated_cardiomyopathy_PES/FEV1_PES_Dilated_cardiomyopathy_all_SNPs_threshold.csv")

########################
## Extension of telomeres PES
########################

## Read in list of genes used to make the pathways in cancer score

Extension_of_telomeres_genes <- read.table("Preparing_geneset_specific_coordinates/Stable_IDs/Telomeres_stable_IDs.txt", header = TRUE)

## Merge pathways in cancer genes with those available in the geuvadis dataset

Merged_Extension_of_telomeres <- merge(Extension_of_telomeres_genes, Available_ID, by = "Stable_ID")

List_Extension_of_telomeres <- as.list(Merged_Extension_of_telomeres)[[1]]

## Define linear model to test the effect of FEV1 Pathways in cancer PES (P < 1) on these genes

PES_Extension_of_telomeres_model <- function(v){
  PES_Extension_of_telomeres_regresssion_result <- lm(glue::glue('{v} ~ SEX + PC1 + PC2 + PC3 + Genome_wide_PGS_0_05_threshold + Extension_of_telomeres_PES'), data = Merged_genotype_mRNA_geuvadis)
  return(summary(PES_Extension_of_telomeres_regresssion_result)) 
}

## Iterate function over list of genes

PES_Extension_of_telomeres <- sapply(List_Extension_of_telomeres, PES_Extension_of_telomeres_model)

## Extract t-statistic and p value for each gene

PES_Extract_Extension_of_telomeres <- apply(PES_Extension_of_telomeres, 2, function(x) return(as.data.frame(x$coefficients)[7,3:4]))

PES_Results_Extension_of_telomeres <- data.frame()
for (i in 1:length(PES_Extract_Extension_of_telomeres)) {
  PES_Results_Extension_of_telomeres <- rbind(PES_Results_Extension_of_telomeres, PES_Extract_Extension_of_telomeres[[i]])
}
rownames(PES_Results_Extension_of_telomeres) <- List_Extension_of_telomeres

## Apply FDR correction and export results

PES_Extension_of_telomeres_FDR <- p.adjust(PES_Results_Extension_of_telomeres$'Pr(>|t|)', method = "fdr")
PES_Extension_of_telomeres_FINAL <- cbind(PES_Results_Extension_of_telomeres, PES_Extension_of_telomeres_FDR)
write.csv(PES_Extension_of_telomeres_FINAL, file = "Extension_of_telomeres_PES/FEV1_PES_Extension_of_telomeres_0_05_SNPs_threshold.csv")

remove(list = ls())
