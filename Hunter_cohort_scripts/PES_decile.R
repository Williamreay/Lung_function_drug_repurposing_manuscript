#####################

## Association between bottom decile PES for class b2 secretin and normalised FVC

## William Reay - May 2020 - github: https://github.com/Williamreay

#####################

setwd("~/Desktop/Pneumonia_cytokine_lung_function/HCS_PES_profiles/")

## Load dependencies

library(dplyr)
library(readxl)
library(MASS)
library(ggplot2)
library(ggpubr)

## Load Phenotype data (normalised residuals for FEV1 and FVC)

HCS_pheno <- read.csv("Spirometry_phenotype_data/All_participants/All_HCS_participants_normalised_spirometry.csv", header = TRUE, na.strings = ' ')

FVC_HCS_PES_combined <- read_excel("FVC_HCS_PES_combined.xlsx")

## Load conversion between Hxxxx IDs and EIDs

Electoral_HCSID_conversion <- read_excel("Electoral_HCSID_conversion.xlsx")

### FVC class b2 secretin PES ####

## Derive EIDs from Hxxxx IDs

EIDs_FVC <- merge(FVC_HCS_PES_combined, Electoral_HCSID_conversion, by = "ID")

## Combine with phenotype data

PES_stats_FVC <- merge(EIDs_FVC, HCS_pheno, by = "electoralId")

## Create a dummy variable for individuals with a bottom decile class b2 secretin PES

PES_stats_FVC$Bottom_decile_b2_secretin_PES <- ntile(PES_stats_FVC$Class_b2_Secretin_PES, 10)

PES_stats_FVC$Bottom_decile_b2_secretin_PES <- ifelse(PES_stats_FVC$Bottom_decile_b2_secretin_PES == 1, 1, 0)

## Test association between normalised FVC and bottom decile PES

Bottom_decile_class_b2_secretin <- glm(Bottom_decile_b2_secretin_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Normalised_FVC_residuals, family = "binomial", data = PES_stats_FVC)

exp(cbind(coef(Bottom_decile_class_b2_secretin), confint(Bottom_decile_class_b2_secretin)))  

## Create a dummy variable for individuals with a bottom decile PGS (P < 0.005)

PES_stats_FVC$Bottom_decile_PGS_0.005 <- ntile(PES_stats_FVC$Genome_wide_0.005_threshold, 10)

PES_stats_FVC$Bottom_decile_PGS_0.005 <- ifelse(PES_stats_FVC$Bottom_decile_PGS_0.005 == 1, 1, 0)

## Test association between normalised FVC and bottom decile PGS

Bottom_decile_PGS <- glm(Bottom_decile_PGS_0.005 ~ PC1 + PC2 + PC3 + PC4 + PC5 + Normalised_FVC_residuals, family = "binomial", data = PES_stats_FVC)
  
exp(cbind(coef(Bottom_decile_PGS), confint(Bottom_decile_PGS)))


