#####################################
## Association between lung function PES and FEV1 or FVC (PES nominally associated with either variable from PRSice analyses)
## William Reay - May 2020
#####################################

setwd("~/Desktop/Pneumonia_cytokine_lung_function/HCS_PES_profiles/")

## Load dependencies
library(readxl)

## Load Phenotype data (normalised residuals for FEV1 and FVC)

HCS_pheno <- read.csv("Spirometry_phenotype_data/All_participants/All_HCS_participants_normalised_spirometry.csv", header = TRUE)

## Load PES, PGS, and PC covariate data

FEV1_HCS_PES_combined <- read_excel("FEV1_HCS_PES_combined.xlsx")
FVC_HCS_PES_combined <- read_excel("FVC_HCS_PES_combined.xlsx")

## Load conversion between Hxxxx IDs and EIDs

Electoral_HCSID_conversion <- read_excel("Electoral_HCSID_conversion.xlsx")

### FEV1 ####

## Derive EIDs from Hxxxx IDs

EIDs_FEV1 <- merge(FEV1_HCS_PES_combined, Electoral_HCSID_conversion, by = "ID")

## Combine with phenotype data

PES_stats_FEV1 <- merge(EIDs_FEV1, HCS_pheno, by = "electoralId")

## Scale the PES and PGS

PES_stats_FEV1$Genome_wide_all_SNPs_threshold <- as.numeric(scale(PES_stats_FEV1$Genome_wide_all_SNPs_threshold))
PES_stats_FEV1$Pathways_in_cancer_PES <- as.numeric(scale(PES_stats_FEV1$Pathways_in_cancer_PES))


## Test association between pathways in cancer PES, adjusted for the PGS at the same threshold

FEV1_pathways_in_cancer_NULL <- glm(Normalised_FEV1_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, family = "gaussian", data = PES_stats_FEV1)

FEV1_pathways_in_cancer_FULL <- glm(Normalised_FEV1_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold + Pathways_in_cancer_PES, family = "gaussian", data = PES_stats_FEV1)

### FVC ####

## Derive EIDs from Hxxxx IDs

EIDs_FVC <- merge(FVC_HCS_PES_combined, Electoral_HCSID_conversion, by = "ID")

## Combine with phenotype data

PES_stats_FVC <- merge(EIDs_FVC, HCS_pheno, by = "electoralId")

## Scale PGS and PES

PES_stats_FVC$Genome_wide_all_SNPs_threshold <- as.numeric(scale(PES_stats_FVC$Genome_wide_all_SNPs_threshold))
PES_stats_FVC$Genome_wide_0.05_threshold <- as.numeric(scale(PES_stats_FVC$Genome_wide_0.05_threshold))
PES_stats_FVC$Genome_wide_0.005_threshold <- as.numeric(scale(PES_stats_FVC$Genome_wide_0.005_threshold))
PES_stats_FVC$Circadian_clock_PES <- as.numeric(scale(PES_stats_FVC$Circadian_clock_PES))
PES_stats_FVC$Class_b2_Secretin_PES <- as.numeric(scale(PES_stats_FVC$Class_b2_Secretin_PES))
PES_stats_FEV1$Pathways_in_cancer_PES <- as.numeric(scale(PES_stats_FVC$Pathways_in_cancer_PES))
PES_stats_FVC$ECM_PES <- as.numeric(scale(PES_stats_FVC$ECM_PES))

## Circadian clock FVC PES

FVC_circadian_clock_NULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.05_threshold, family = "gaussian", data = PES_stats_FVC)

FVC_circadian_clock_FULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.05_threshold + Circadian_clock_PES, family = "gaussian", data = PES_stats_FVC)

## Pathways in cancer FVC PES

FVC_pathways_in_cancer_NULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, family = "gaussian", data = PES_stats_FVC)

FVC_pathways_in_cancer_FULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold + Pathways_in_cancer_PES, family = "gaussian", data = PES_stats_FVC)

## Class b2 secretin FVC PES

FVC_secretin_NULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.005_threshold, family = "gaussian", data = PES_stats_FVC)

FVC_secretin_FULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_0.005_threshold + Class_b2_Secretin_PES, family = "gaussian", data = PES_stats_FVC)

# F test of residual deviance

FVC_secretin_dev <- anova(FVC_secretin_NULL, FVC_secretin_FULL, test = "F")

## NABA matrisome (ECM PES)

FVC_ECM_NULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold, family = "gaussian", data = PES_stats_FVC)

FVC_ECM_FULL <- glm(Normalised_FVC_residuals ~ PC1 + PC2 + PC3 + PC4 + PC5 + Genome_wide_all_SNPs_threshold + ECM_PES, family = "gaussian", data = PES_stats_FVC)


