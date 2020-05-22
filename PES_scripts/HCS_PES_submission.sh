#!/bin/bash

cd /home/control/data/users/william/Cytokine_pneumonia_lung_function/Lung_function_PES

while read GWAS PHENOTYPE PATHWAY THRESHOLD;
do
python3 ../../PES_2020/PES_generation.py \
--pathway_gwasfile HCS_PES/${GWAS} \
--bfile ~/data/users/william/HCS-HRC/plink/Clean_HCS_HRC \
--PRsice2_binary ../../PES_2020/PRSice_linux/PRSice_linux  \
--phenotype $PHENOTYPE \
--pathway $PATHWAY \
--A1 Coded --A2 Non_coded --CHR Chromosome --BP Position_b37 --SNP SNP \
--p_threshold $THRESHOLD \
--statistic beta --P P;
done < HCS_PES_pathways.tab.txt
