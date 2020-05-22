#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#Cigs per day - GTEx v7 Lung weights

for chr in $(seq 1 22);
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats ~/data/users/william/Cytokine_pneumonia_lung_function/Smoking_TWAS/Munged/Cigarettes_smoked_per_day_Munged.sumstats.gz \
--weights  ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/Lung.P01.pos \
--weights_dir ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr ${chr} \
--out ~/data/users/william/Cytokine_pneumonia_lung_function/Smoking_TWAS/Lung_results/Cigs_per_day_lung/Cigs_per_day_Lung_chr${chr}_FUSION
done
