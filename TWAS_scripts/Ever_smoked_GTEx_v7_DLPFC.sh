#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#Ever vs never smoked - CMC DLPFC

for chr in $(seq 1 22);
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats ~/data/users/william/Cytokine_pneumonia_lung_function/Smoking_TWAS/Munged/Ever_smoked_munged.sumstats.gz \
--weights  ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/CMC.BRAIN.RNASEQ.pos \
--weights_dir ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr ${chr} \
--out ~/data/users/william/Cytokine_pneumonia_lung_function/Smoking_TWAS/DLPFC_results/Ever_smoked_DLPFC/Ever_smoked_DLPFC_chr${chr}_FUSION
done
