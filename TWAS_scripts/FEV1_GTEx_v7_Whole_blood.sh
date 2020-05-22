#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#FEV1 - GTEx v7 Whole_blood weights

for chr in $(seq 1 22);
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/Munged_TWAS_input/Munged_FEV1_meta.sumstats.gz \
--weights  ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/Whole_Blood.P01.pos \
--weights_dir ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr ${chr} \
--out ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/GTEx_v7_Whole_blood/FEV1_results/FEV1_Whole_blood_chr${chr}_FUSION
done
