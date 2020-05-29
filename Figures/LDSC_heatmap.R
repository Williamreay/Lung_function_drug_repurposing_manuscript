##################################
## Spirometry LDSC heatmap
## William Reay - June 2020
##################################

set.seed(1234)

## Load dependencies

library(readxl)
library(ComplexHeatmap)
library(circlize)

Heatmap_input <- read_excel("Desktop/Pneumonia_cytokine_lung_function/Lung_function_LDhub/Heatmap_input.xlsx")

## Contruct matrix

LDSC_mat <- as.matrix(Heatmap_input[, -1])

## Define the row names.

rownames(LDSC_mat) = c("Age of first birth", "Anorexia Nervosa", "Asthma", "BMI", "Cigarettes per day",
                       "College completion", "Extreme BMI", "Extreme waist:hip", "Fasting glucose", "Fasting insulin",
                       "Father's age at death", "HDL cholesterol", "Hip circumference", "HOMA-B", "HOMA-IR",
                       "Leptin (adj BMI)", "Leptin", "Lung cancer", "Obesity class 1", "Obesity class 2", "Obesity class 3",
                       "Overweight", "Sitting height ratio", "Triglycerides", "Type 2 Diabetes", "Urate",
                       "Waist circumference", "Waist:hip", "Years of Schooling")

## Transpose matrix

LDSC_transposed <- t(LDSC_mat)


## Construct heatmap

ht = Heatmap(LDSC_transposed, name = "rg", rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = TRUE, show_row_dend = FALSE, clustering_distance_rows = "pearson",
        row_dend_width = unit(1, "cm"))

