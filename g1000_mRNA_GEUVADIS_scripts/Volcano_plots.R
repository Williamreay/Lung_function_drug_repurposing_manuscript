#####################

## Volcano plots for mRNA/PES correlations

## William Reay - May 2020 - github: https://github.com/Williamreay

#####################

library(ggplot2)
library(calibrate)
library(readr)

setwd("~/Desktop/Pneumonia_cytokine_lung_function/mRNA_PES_correlation")

## Class b2 secretin

FVC_PES_Class_b2_secretin_0_005_SNPs_threshold <- read_csv("Class_b2_secretin_PES/FVC_PES_Class_b2_secretin_0_005_SNPs_threshold.csv")

FVC_PES_Class_b2_secretin_0_005_SNPs_threshold$threshold <- factor(ifelse(FVC_PES_Class_b2_secretin_0_005_SNPs_threshold$PES_Class_b2_secretin_FDR < 0.1,1,ifelse(FVC_PES_Class_b2_secretin_0_005_SNPs_threshold $PES_Class_b2_secretin_FDR < 0.05,-1, 0)))


Plot_geuvadis_secretin <- ggplot(data = FVC_PES_Class_b2_secretin_0_005_SNPs_threshold,
                                  aes(x=t_value, y = -log10(P),
                                  colour = threshold, label = Sig_gene_IDs)) +
                                  scale_color_manual(name="FDR", values = c("black","#3333FF")) +
                                  geom_point(alpha = 0.8, size = 1.75) +
                                  labs(x = expression("t value"), y = expression(paste("-log"[10], "P-value"))) +
                                  theme_bw() +
                                  geom_text(aes(label=ifelse(PES_Class_b2_secretin_FDR < 0.05, as.character(Sig_gene_IDs), '')), hjust=-0.1, vjust=-0.5) +
                                  theme(legend.position ="none") +
                                  theme(legend.position ="none") +
                                  theme(axis.title = element_text(face="bold", size=12)) +
                                  ggtitle("FVC:Class B/2 secretin family receptors") +
                                  xlim(c(-4,4)) +
                                  geom_hline(yintercept = 1.3, linetype ="longdash") +
                                  ylim(c(0, 5))


## Circadian clock

FVC_PES_Circadian_clock_0_05_SNPs_threshold <- read_csv("Circadian_clock_PES/FVC_PES_Circadian_clock_0_05_SNPs_threshold.csv")

FVC_PES_Circadian_clock_0_05_SNPs_threshold$threshold <- factor(ifelse(FVC_PES_Circadian_clock_0_05_SNPs_threshold$PES_Circadian_clock_FDR < 0.1,1,ifelse(FVC_PES_Circadian_clock_0_05_SNPs_threshold$PES_Circadian_clock_FDR < 0.05,-1, 0)))


Plot_geuvadis_circadian <- ggplot(data = FVC_PES_Circadian_clock_0_05_SNPs_threshold,
                                 aes(x=t_value, y = -log10(P),
                                 colour = threshold, label = Sig_gene_IDs)) +
                                    geom_point(alpha = 0.8, size = 1.75) +
                                    scale_color_manual(name="FDR", values = c("black","red")) +
                                    labs(x = expression("t value"), y = expression(paste("-log"[10], "P-value"))) +
                                    geom_text(aes(label=ifelse(PES_Circadian_clock_FDR < 0.1, as.character(Sig_gene_IDs), '')), hjust=1, vjust=-0.5) +
                                    theme(legend.position ="none") +
                                    theme_bw() +
                                    theme(legend.position ="none") +
                                    theme(axis.title = element_text(face="bold", size=12)) +
                                    ggtitle("FVC:Circadian clock") +
                                    xlim(c(-4,4)) +
                                    geom_hline(yintercept = 1.3, linetype ="longdash") +
                                    ylim(c(0, 5))

## Pathways in cancer

FVC_PES_pathways_in_cancer_all_SNPs_threshold <- read_csv("Pathways_in_cancer_PES/FVC_PES_pathways_in_cancer_all_SNPs_threshold.csv")

FVC_PES_pathways_in_cancer_all_SNPs_threshold$threshold <- factor(ifelse(FVC_PES_pathways_in_cancer_all_SNPs_threshold$PES_pathways_in_cancer_FDR < 0.1,1,ifelse(FVC_PES_pathways_in_cancer_all_SNPs_threshold$PES_pathways_in_cancer_FDR < 0.05,-1, 0)))


Plot_geuvadis_cancer <- ggplot(data = FVC_PES_pathways_in_cancer_all_SNPs_threshold,
                                  aes(x=t_value, y = -log10(P),
                                    colour = threshold, label = Sig_gene_IDs)) +
                                    geom_point(alpha = 0.8, size = 1.75) +
                                    scale_color_manual(name="FDR", values = c("black","red")) +
                                    labs(x = expression("t value"), y = expression(paste("-log"[10], "P-value"))) +
                                    geom_text(aes(label=ifelse(PES_pathways_in_cancer_FDR < 0.1, as.character(Sig_gene_IDs), '')), hjust=1, vjust=-0.5) +
                                    theme(legend.position ="none") +
                                    theme_bw() +
                                    theme(legend.position ="none") +
                                    theme(axis.title = element_text(face="bold", size=12)) +
                                    ggtitle("FVC:Pathways in cancer") +
                                    xlim(c(-4,4)) +
                                    geom_hline(yintercept = 1.3, linetype ="longdash") +
                                    ylim(c(0, 5))
