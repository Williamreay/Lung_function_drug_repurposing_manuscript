#############################
## Effect of fasting glucose on FVC
## William Reay 2020
#############################

setwd("~/Desktop/Pneumonia_cytokine_lung_function/Urate_glucose_MR/")

## Load dependencies
library(dplyr)
library(TwoSampleMR)
library(MRPRESSO)

## Import IVs for Fasting glucose (ln mmol/L) and clump for LD
Glucose_IVs <- read_exposure_data("Fasting_glucose_Scott_et_al.tab.txt", clump = TRUE, sep = "\t")

## Load R data object with FVC GWAS summary statistics in MR base format
load("../Cytokine_MR/FVC.RData")

## Extract instrument variables from FVC GWAS
Glucose_FVC_outcome <- format_data(
  dat = FVC_GWAS,
  type = "outcome",
  snps = Glucose_IVs$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf")

## Harmonise data and exclude palindromic SNPs
Glucose_FVC_harmonised <- harmonise_data(exposure_dat = Glucose_IVs, 
                                                   outcome_dat = Glucose_FVC_outcome, 
                                                   action = 3)

## Perform MR (TwoSampleMR included tests - IVW-MRE, Weighted Median, MR-Egger)
Glucose_to_FVC <- mr(Glucose_FVC_harmonised ,method_list=c("mr_ivw_mre", 
                                                           "mr_ivw_fe",
                                                             "mr_egger_regression",
                                                              "mr_weighted_median"))
## Get beta 95% CI
generate_odds_ratios(Glucose_to_FVC)

## Convert SE to SD and perform MR-PRESSO

Filtered_glucose_FVC <- Glucose_FVC_harmonised %>% filter(mr_keep == TRUE)

MR_PRESSO_input_glucose_FVC <- Filtered_glucose_FVC[,c("SNP", "beta.outcome", "beta.exposure", "se.outcome", "se.exposure")]

MR_PRESSO_glucose_FVC <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
                                          DISTORTIONtest = TRUE, data = MR_PRESSO_input_glucose_FVC,
                                          NbDistribution = 10000)

## Heterogeneity test via Cochran's Q
Glucose_FVC_het <- mr_heterogeneity(Glucose_FVC_harmonised)

## Test if Egger intercept is significantly different from zero
Glucose_FVC_Egger_intercept <- mr_pleiotropy_test(Glucose_FVC_harmonised)

## Leave one out analysis to test for IVs which disproprtionately contribute to the causal estimate
Glucose_FVC_LOO <- mr_leaveoneout(Glucose_FVC_harmonised)

Glucose_FVC_single_SNP <- mr_singlesnp(Glucose_FVC_harmonised)

Glucose_FVC_LOO_plot <- mr_leaveoneout_plot(Glucose_FVC_LOO)

## Scatter plot Glucose ==> FVC
Glucose_FVC_sp <- mr_scatter_plot(Glucose_to_FVC, Glucose_FVC_harmonised)
