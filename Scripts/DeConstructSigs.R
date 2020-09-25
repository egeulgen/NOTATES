##################################################
## Project: NOTATES
## Script purpose: Script for identifying COSMIC
## mutational signatures (SBS v3) in the sample
## using deconstructSigs
## Date: Sep 23, 2020
## Author: Ege Ulgen
##################################################

# load and install required libraries -------------------------------------
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

if(!suppressPackageStartupMessages(require(deconstructSigs))) {
  devtools::install_github("raerose01/deconstructSigs")
  suppressPackageStartupMessages(library(deconstructSigs))
}

arguments <- commandArgs(trailingOnly = T)
script_dir <- arguments[1]
sample_id <- arguments[2]

# load somatic SNVs -------------------------------------------------------
somatic_SNVs <- read.delim("./Oncotator/annotated.sSNVs.tsv", comment.char="#")

# Subsetting for (MuTect's default) HQ filters
somatic_SNVs <- subset(somatic_SNVs, 
                       alt_allele_seen=="True")
# filter for VAF > 0.05
somatic_SNVs <- somatic_SNVs[somatic_SNVs$tumor_f > 0.05, ]


somatic_SNVs$Sample <- sample_id
somatic_SNVs <- somatic_SNVs[,c("Sample","Chromosome", "Start_position","Reference_Allele", "Tumor_Seq_Allele2")]
somatic_SNVs$Chromosome <- paste0("chr", somatic_SNVs$Chromosome)

# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = somatic_SNVs, 
                                sample.id = "Sample", 
                                chr = "Chromosome", 
                                pos = "Start_position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2")

# Determine the signatures contributing to the sample
signatures <- whichSignatures(tumor.ref = sigs.input,
                              signatures.ref = signatures.exome.cosmic.v3.may2019,
                              sample.id = sample_id,
                              contexts.needed = TRUE,
                              tri.counts.method = 'exome')

# Plots of output
pdf("mut_signatures.pdf")
plotSignatures(signatures, sub = '')
dev.off()

pdf("mut_signatures_pie.pdf")
makePie(signatures, v3 = TRUE)
dev.off()

# Add details (manually extracted from COSMIC)
nms <- names(signatures$weights)[signatures$weights!=0]
cosmic_sigs <- read.csv(paste0(script_dir,"/cosmic_sbs_details.csv"), stringsAsFactors = F)

sig_df <- t(round(as.matrix(signatures$weights[nms]),2))
rownames(sig_df) <- gsub("\\.", " ",rownames(sig_df))
colnames(sig_df) <- "freq"
sig_df <- cbind(sig_df,cosmic_sigs[match(rownames(sig_df), cosmic_sigs$Signature), 2:4])

write.csv(sig_df, "mut_signatures.csv")
