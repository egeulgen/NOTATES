# load and install required libraries -------------------------------------
if(!require(BSgenome.Hsapiens.UCSC.hg19))
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
}
if(!require(deconstructSigs))
{
  if(!require(devtools))
    install.packages("devtools")
  devtools::install_github("raerose01/deconstructSigs")
}
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg19)

arguments <- commandArgs(trailingOnly = T)
script_dir <- arguments[1]
sample_id <- arguments[2]

# load somatic SNVs -------------------------------------------------------
somatic_SNVs <- read.delim("./Oncotator/annotated.sSNVs.tsv", stringsAsFactors=F, comment.char="#")

# Subsetting for (MuTect's default) HQ filters
somatic_SNVs <- subset(somatic_SNVs, 
                       alt_allele_seen=="True")
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
makePie(signatures)
dev.off()

nms <- names(signatures$weights)[signatures$weights!=0]

cosmic_sigs <- read.csv(paste0(script_dir,"/cosmic_sbs_details.csv"), stringsAsFactors = F)

sig_df <- t(round(as.matrix(signatures$weights[nms]),2))
rownames(sig_df) <- gsub("\\.", " ",rownames(sig_df))
colnames(sig_df) <- "freq"
sig_df <- cbind(sig_df,cosmic_sigs[match(rownames(sig_df), cosmic_sigs$Signature), 2:4])

write.csv(sig_df, "mut_signatures.csv")
