##################################################
## Project: NOTATES
## Script purpose: Prepare MAF file for MSIpred
## Analysis
## Date: Sep 22, 2020
## Author: Ege Ulgen
##################################################

args <- commandArgs(trailingOnly=TRUE)

dir.create("MSIpred")
funcotator_maf <- read.delim("Funcotator/annotated_somatic.maf", comment.char="#")
funcotator_maf$Tumor_Sample_Barcode <- args[1]

nec_cols <- c("Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode")
funcotator_maf <- funcotator_maf[, nec_cols]

chr.list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
funcotator_maf <- funcotator_maf[funcotator_maf$Chromosome %in% chr.list, ]

write.table(funcotator_maf, "MSIpred/somatic_maf.maf", sep = "\t", quote = FALSE, row.names = FALSE) 
