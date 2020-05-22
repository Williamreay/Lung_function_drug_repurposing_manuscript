########
## Disentangling correlation from causation in the putatitve relationship between significantly genetically correlated traits and FEV1/FVC
## 2020 - William Reay
########

setwd("~/Desktop/Pneumonia_cytokine_lung_function/LCV")

##Load dependencies
suppressMessages(library(optparse))

##Specify command line inputs
option_list = list(
  make_option("--phenotype_name", action="store", default=NA, type='character',
              help="The name of the trait to be analysed with FEV1:FVC [required]"),
  make_option("--sumstats_file", action="store", default=NA, type='character',
              help="File name for summary stats for trait to be analysed, gzip compression required")
)
opt = parse_args(OptionParser(option_list=option_list))

#Load in munged summary statistics for FEV1/FVC and trait to be analysed

cat("#########################")
cat("\n")
cat("Loading summary statistics for FEV1/FVC and", paste(opt$phenotype_name))
cat("\n")
cat("#########################")
cat("\n")

#Load trait 1 data 
opt$phenotype_name_df <- read.table(file = paste("Munged/",opt$sumstats_file,sep=""),
                                    header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))

opt$phenotype_name_df <- na.omit(opt$phenotype_name_df)

#Load trait 2 data 
FVC_df <- read.table("../LDSC/LD_hub/Munged_FEV1_FVC_meta.sumstats.gz" ,header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))
FVC_df <- na.omit(FVC_df)

cat("\n")
cat("#########################")
cat("\n")
cat("GWAS data loaded")
cat("\n")
cat("#########################")
cat("\n")

cat("\n")
cat("#########################")
cat("\n")
cat("Loading LD scores")
cat("\n")
cat("#########################")
cat("\n")

#Load LD scores
LD_scores <- read.table("unannotated_LDscore.l2.ldsc",header=TRUE,sep="\t",stringsAsFactors=FALSE)

cat("\n")
cat("#########################")
cat("\n")
cat("Merging data frames")
cat("\n")
cat("#########################")
cat("\n")


#Merge such that SNPs are annotated with LD scores
Merged_df <- merge(LD_scores,opt$phenotype_name_df,by="SNP")
Annotated_df_combined <- merge(Merged_df,FVC_df,by="SNP")

#Sort by position 
Sorted_df <- Annotated_df_combined[order(Annotated_df_combined[,"CHR"],Annotated_df_combined[,"BP"]),]

Sorted_df$L2 <- as.numeric(Sorted_df$L2)

#Check if any mismatches
mismatch = which(Sorted_df$A1.x!=Sorted_df$A1.y,arr.ind=TRUE)
Sorted_df[mismatch,]$Z.y = Sorted_df[mismatch,]$Z.y*-1
Sorted_df[mismatch,]$A1.y = Sorted_df[mismatch,]$A1.x
Sorted_df[mismatch,]$A2.y = Sorted_df[mismatch,]$A2.x


cat("\n")
cat("#########################")
cat("\n")
cat("Constructing LCV model, output will be sent to text file",  paste(opt$phenotype_name, '_LCV_FVC.txt', sep = ""))
cat("\n")
cat("#########################")
cat("\n")


#Run LCV-need to setwd to directory containing LCV package
source("RunLCV.R")

file = paste(opt$phenotype_name, '_LCV_FEV1_FVC.txt', sep = "")

sink(file, append = TRUE)

LCV <- RunLCV(Sorted_df$L2,Sorted_df$Z.x, Sorted_df$Z.y)
sprintf("Estimated posterior gcp=%.2f(%.2f), pvalue=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, LCV$pval.gcpzero.2tailed, LCV$rho.est, LCV$rho.err)

#Sink output to text file

LCV

sink()

#Clear environment after run
rm(list = ls())