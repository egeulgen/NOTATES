# load and install required libraries -------------------------------------
if(!require(deconstructSigs))
{
  if(!require(devtools))
    install.packages("devtools")
  install_github("raerose01/deconstructSigs")
}
library(deconstructSigs)
if(!require(BSgenome.Hsapiens.UCSC.hg19))
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)

arguments <- commandArgs(trailingOnly = T)
script_dir <- arguments[1]

# load somatic SNVs -------------------------------------------------------
somatic_SNVs <- read.delim("./Oncotator/annotated.sSNVs.tsv", stringsAsFactors=F, comment.char="#")

# Subsetting for (MuTect's default) HQ filters
somatic_SNVs <- subset(somatic_SNVs, 
                       alt_allele_seen=="True" & 
                         alt_allele_in_normal=="PASS" & 
                         clustered_events=="PASS" &
                         germline_risk == "PASS" & 
                         homologous_mapping_event=="PASS" &
                         multi_event_alt_allele_in_normal=="PASS" &
                         panel_of_normals=="PASS" & 
                         str_contraction=="PASS" &
                         t_lod_fstar=="PASS" &
                         triallelic_site=="PASS" &
                         short_tandem_repeat_membership=="False")
somatic_SNVs$tumor <- "tumor"
somatic_SNVs <- somatic_SNVs[,c("tumor","Chromosome", "Start_position","Reference_Allele", "Tumor_Seq_Allele2")]

# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = somatic_SNVs, 
                                sample.id = "tumor", 
                                chr = "Chromosome", 
                                pos = "Start_position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2")

# Determine the signatures contributing to the sample
signatures <- whichSignatures(tumor.ref = sigs.input,
                              signatures.ref = signatures.cosmic,
                              sample.id = "tumor",
                              contexts.needed = TRUE,
                              tri.counts.method = 'default')

# Plots of output
pdf("mut_signatures.pdf")
plotSignatures(signatures, sub = '')
dev.off()
pdf("mut_signatures_pie.pdf")
nms <- names(signatures$weights)[signatures$weights!=0]
pie(as.matrix(signatures$weights[nms]))
dev.off()

cosmic_sigs <- read.csv(paste0(script_dir,"/COSMIC_signatures.csv"), stringsAsFactors = F)

sig_df <- t(round(as.matrix(signatures$weights[nms]),2))
rownames(sig_df) <- gsub("\\.", " ",rownames(sig_df))
colnames(sig_df) <- "freq"
sig_df <- cbind(sig_df,cosmic_sigs[match(rownames(sig_df),cosmic_sigs$Signature),2:5])

write.csv(sig_df, "mut_signatures.csv")
