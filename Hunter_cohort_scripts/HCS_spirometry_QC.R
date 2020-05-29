#####################

## Cleaning and normalisation of the Hunter Community Study cohort spirometry data

## William Reay - May 2020 - github: https://github.com/Williamreay

#####################

## Load dependencies

library(dplyr)
library(RNOmni)
library(easyGgplot2)
library(MASS)

## Import HCS phenotype data

Raw_HCS <- read.csv("~/Desktop/Hunter_cohort/Phenotype_data/mcv1.csv", header = TRUE, na.strings = ' ')


## Remove individuals with missing maximum FEV1 (fev1w1) and FVC (fvcw1)

Raw_spirometry_HCS <- Raw_HCS %>% filter(fev1w1 != 'NA' & fvcw1 != 'NA')


## Plot histogram of untransformed values for both spirometry indices 

Histo_FEV1 <- ggplot2.histogram(data = Raw_spirometry_HCS, xName="fev1w1",
                                fill = "white", color = "black", addMeanLine = TRUE,
                                meanLineType = "dashed", meanLineColor = "red",
                                addDensityCurve = TRUE, densityFill = '#FF6666',
                                xtitle = "Maximum FEV1", bins = 50)

Histo_FVC <- ggplot2.histogram(data = Raw_spirometry_HCS, xName="fvcw1",
                                fill = "white", color = "black", addMeanLine = TRUE,
                                meanLineType = "dashed", meanLineColor = "red",
                                addDensityCurve = TRUE, densityFill = '#0066FF',
                                xtitle = "Maximum FVC", bins = 50)

## Remove individuals who are missing any of the following values:
## a) Height, b) Sex, c) Age, d) Asthma status, e) Bronchitis/emphysema status

Spirometry_regression_input_HCS <- Raw_spirometry_HCS %>% filter(SmokeEverw1 != 'NA' & Heightw1 != 'NA' & Bronchitis_Emphysemaw1 != 'NA' & Asthmaw1 != 'NA' & sex != 'NA' & bage != 'NA')

## Linear model to test the effect of the above covariates on FEV1 and FVC respectively (test for quadratic effects of age and height)

FEV1_regression <- lm(fev1w1 ~ sex + bage + I(bage^2) + Heightw1 + I(Heightw1^2) + SmokeEverw1 + Asthmaw1 + Bronchitis_Emphysemaw1, data = Spirometry_regression_input_HCS)

FVC_regression <- lm(fvcw1 ~ sex + bage + I(bage^2) + Heightw1 + I(Heightw1^2) + SmokeEverw1 + Asthmaw1 + Bronchitis_Emphysemaw1, data = Spirometry_regression_input_HCS)

## Extract resdiuals from both models

FEV1_residuals_RAW <- as.data.frame(FEV1_regression$residuals)
FVC_residuals_RAW <- as.data.frame(FVC_regression$residuals)

colnames(FEV1_residuals_RAW)[1] <- "FEV1_residuals"
colnames(FVC_residuals_RAW)[1] <- "FVC_residuals"

## Inverse-rank normalisation of the residuals  - use default offset (Blom transformation)

FEV1_residuals_NORMALISED <- as.data.frame(rankNorm(FEV1_residuals_RAW$FEV1_residuals), k=3/8)

FVC_residuals_NORMALISED <- as.data.frame(rankNorm(FVC_residuals_RAW$FVC_residuals), k=3/8)

colnames(FEV1_residuals_NORMALISED)[1] <- "Normalised_FEV1_residuals"
colnames(FVC_residuals_NORMALISED)[1] <- "Normalised_FVC_residuals"

## Plot histogram of transformed values for both spirometry indices 

Histo_FEV1_transformed <- ggplot2.histogram(data = FEV1_residuals_NORMALISED, xName="Normalised_FEV1_residuals",
                                fill = "white", color = "black", addMeanLine = TRUE,
                                meanLineType = "dashed", meanLineColor = "red",
                                addDensityCurve = TRUE, densityFill = '#FF6666',
                                xtitle = "Normalised maximum FEV1", bins = 50)

Histo_FVC_transformed <- ggplot2.histogram(data = FVC_residuals_NORMALISED, xName="Normalised_FVC_residuals",
                               fill = "white", color = "black", addMeanLine = TRUE,
                               meanLineType = "dashed", meanLineColor = "red",
                               addDensityCurve = TRUE, densityFill = '#0066FF',
                               xtitle = "Normalised maximum FVC", bins = 50)

## Merge normalised residuals with regression input df

Spirometry_regression_input_HCS$Normalised_FEV1_residuals <- FEV1_residuals_NORMALISED$Normalised_FEV1_residuals
Spirometry_regression_input_HCS$Normalised_FVC_residuals <- FVC_residuals_NORMALISED$Normalised_FVC_residuals

## Write output to .csv

write.csv(Spirometry_regression_input_HCS, file = "~/Desktop/Pneumonia_cytokine_lung_function/HCS_PES_profiles/Spirometry_phenotype_data/All_participants/All_HCS_participants_normalised_spirometry.csv",
          row.names = FALSE, quote = FALSE)

## Extract individuals who were both a) non-smokers and b) have no respiratory illnesses (asthma, emphysema, or bronchitis)

Smokers_respiratory_removed_regression_input <- Raw_spirometry_HCS %>% filter(SmokeEverw1 == 0 & Bronchitis_Emphysemaw1 == 0 & Asthmaw1 == 0 & SmokeEverw1 != 'NA' & Heightw1 != 'NA' & Bronchitis_Emphysemaw1 != 'NA' & Asthmaw1 != 'NA' & sex != 'NA' & bage != 'NA')

## Linear model to test the effect of the above covariates on FEV1 and FVC respectively (test for quadratic effects of age and height)

Non_smokers_FEV1_regression <- lm(fev1w1 ~ sex + bage + I(bage^2) + Heightw1 + I(Heightw1^2), data = Smokers_respiratory_removed_regression_input)

Non_smokers_FVC_regression <- lm(fvcw1 ~ sex + bage + I(bage^2) + Heightw1 + I(Heightw1^2), data = Smokers_respiratory_removed_regression_input)

## Extract resdiuals from both models

Non_smokers_FEV1_residuals_RAW <- as.data.frame(Non_smokers_FEV1_regression$residuals)
Non_smokers_FVC_residuals_RAW <- as.data.frame(Non_smokers_FVC_regression$residuals)

colnames(Non_smokers_FEV1_residuals_RAW)[1] <- "FEV1_residuals"
colnames(Non_smokers_FVC_residuals_RAW)[1] <- "FVC_residuals"

## Inverse-rank normalisation of the residuals  - use default offset (Blom transformation)

Non_smokers_FEV1_residuals_NORMALISED <- as.data.frame(rankNorm(Non_smokers_FEV1_residuals_RAW$FEV1_residuals), k=3/8)

Non_smokers_FVC_residuals_NORMALISED <- as.data.frame(rankNorm(Non_smokers_FVC_residuals_RAW$FVC_residuals), k=3/8)

colnames(Non_smokers_FEV1_residuals_NORMALISED)[1] <- "Normalised_FEV1_residuals"
colnames(Non_smokers_FVC_residuals_NORMALISED)[1] <- "Normalised_FVC_residuals"

## Merge normalised residuals with regression input df

Smokers_respiratory_removed_regression_input$Normalised_FEV1_residuals <- Non_smokers_FEV1_residuals_NORMALISED$Normalised_FEV1_residuals
Smokers_respiratory_removed_regression_input$Normalised_FVC_residuals <- Non_smokers_FVC_residuals_NORMALISED$Normalised_FVC_residuals

## Write output to .csv

write.csv(Smokers_respiratory_removed_regression_input, file = "~/Desktop/Pneumonia_cytokine_lung_function/HCS_PES_profiles/Spirometry_phenotype_data/Smokers_respiratory_illness_excluded/Respiratory_illness_smokers_excluded_HCS_participants_normalised_spirometry.csv",
          row.names = FALSE, quote = FALSE)

## Plot histograms before and after normalisation for this cohort

Non_smokers_Histo_FEV1_raw <- ggplot2.histogram(data = Smokers_respiratory_removed_regression_input, xName="fev1w1",
                                            fill = "white", color = "black", addMeanLine = TRUE,
                                            meanLineType = "dashed", meanLineColor = "red",
                                            addDensityCurve = TRUE, densityFill = '#FF6666',
                                            xtitle = "Maximum FEV1", bins = 50)

Non_smokers_Histo_FVC_raw <- ggplot2.histogram(data = Smokers_respiratory_removed_regression_input, xName="fvcw1",
                                           fill = "white", color = "black", addMeanLine = TRUE,
                                           meanLineType = "dashed", meanLineColor = "red",
                                           addDensityCurve = TRUE, densityFill = '#0066FF',
                                           xtitle = "Maximum FVC", bins = 50)

Non_smokers_Histo_FEV1_normalised <- ggplot2.histogram(data = Smokers_respiratory_removed_regression_input, xName="Normalised_FEV1_residuals",
                                                fill = "white", color = "black", addMeanLine = TRUE,
                                                meanLineType = "dashed", meanLineColor = "red",
                                                addDensityCurve = TRUE, densityFill = '#FF6666',
                                                xtitle = "Normalised maximum FEV1", bins = 50)

Non_smokers_Histo_FVC_normalised <- ggplot2.histogram(data = Smokers_respiratory_removed_regression_input, xName="Normalised_FVC_residuals",
                                               fill = "white", color = "black", addMeanLine = TRUE,
                                               meanLineType = "dashed", meanLineColor = "red",
                                               addDensityCurve = TRUE, densityFill = '#0066FF',
                                               xtitle = "Normalised maximum FVC", bins = 50)
