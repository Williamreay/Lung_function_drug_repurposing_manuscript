########
## Glucose - Cigarettes per day LCV
## 09/01/20
########

setwd("~/Desktop/Pneumonia_cytokine_lung_function/LCV/")

#Load in munged summary statistics for Anorexia and fasting insulin (BMI adjusted) - note that the Manning et al summary statistics are used here instead of Scott et al, as insufficient SNPs were available in the MetaboChip summary statistics from that publication

FG="Munged/Munged_fasting_glucose.sumstats.gz"

Cigs="Munged/Munged_hapmap_merged_cigs_per_day.sumstats.gz"

#Load trait 1 data and calculate Zs
FG_df = na.omit(read.table(gzfile(FG),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load trait 2 data and calculate Zs
Cigs_df = na.omit(read.table(gzfile(Cigs),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load LD scores
LD_scores=read.table("unannotated_LDscore.l2.ldsc",header=TRUE,sep="\t",stringsAsFactors=FALSE)

#Merge such that SNPs are annotated with LD scores
Merged_df = merge(LD_scores,FG_df,by="SNP")
Annotated_df_FG_Cigs = merge(Merged_df,Cigs_df,by="SNP")

#Sort by position 
Sorted_df = Annotated_df_FG_Cigs[order(Annotated_df_FG_Cigs[,"CHR"],Annotated_df_FG_Cigs[,"BP"]),]

Sorted_df$L2 = as.numeric(Sorted_df$L2)

#Check if any mismatches
mismatch = which(Sorted_df$A1.x!=Sorted_df$A1.y,arr.ind=TRUE)
Sorted_df[mismatch,]$Z.y = Sorted_df[mismatch,]$Z.y*-1
Sorted_df[mismatch,]$A1.y = Sorted_df[mismatch,]$A1.x
Sorted_df[mismatch,]$A2.y = Sorted_df[mismatch,]$A2.x

#Run LCV-need to setwd to directory containing LCV package
source("RunLCV.R")

LCV = RunLCV(Sorted_df$L2,Sorted_df$Z.x, Sorted_df$Z.y)
sprintf("Estimated posterior gcp=%.2f(%.2f), pvalue=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, LCV$pval.gcpzero.2tailed, LCV$rho.est, LCV$rho.err)
