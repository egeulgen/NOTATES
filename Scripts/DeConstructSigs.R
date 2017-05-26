# load and install required libraries -------------------------------------
if(!require(deconstructSigs))
{
  if(!require(devtools))
    install.packages("devtools")
  if(!require(BSgenome.Hsapiens.UCSC.hg19))
  {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg19")
  }
  install_github("raerose01/deconstructSigs")
}
library(BSgenome.Hsapiens.UCSC.hg19)

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

# Plot output
pdf("mut_signatures.pdf")
plotSignatures(signatures, sub = '')
nms <- names(signatures$weights)[signatures$weights!=0]
barplot(as.matrix(signatures$weights[nms]), ylim = c(0,1))
dev.off()

write.csv(round(as.matrix(signatures$weights[nms]),2), "mut_signatures.csv", row.names = F)
