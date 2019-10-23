##################################################
## Project: NOTATES
## Script purpose: Prepare MAF file for MSIpred
## Analysis
## Date: Oct 22, 2019
## Author: Ege Ulgen
##################################################
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)

dir.create("MSIpred")
oncotator_maf <- read.delim("Oncotator/annotated.sSNVs.tsv", comment.char="#")
oncotator_maf$Chromosome <- paste0("chr", oncotator_maf$Chromosome)
colnames(oncotator_maf) <- sub("position", "Position", colnames(oncotator_maf))
oncotator_maf$Tumor_Sample_Barcode <- args[1]
write.table(oncotator_maf, "MSIpred/somatic_maf.maf", sep = "\t", quote = FALSE, row.names = FALSE) 
