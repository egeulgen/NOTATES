############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                    ExomeCNV Analysis                     #
#             GATK pileup files to baf files               #
#                   Ege Ulgen, Oct 2016                    #
############################################################

setwd("./ExomeCNV/baf")

# read VCF files -----------------------------------------------
cnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample")
normal_VCF <- read.delim("./normal_HQ_SNPs.vcf", header = F, col.names = cnames, comment.char = "#", stringsAsFactors = F)
tumor_VCF <- read.delim("./tumor_HQ_SNPs.vcf",  header = F, col.names = cnames, comment.char = "#", stringsAsFactors = F)

# parse and add coverages -------------------------------------------------
normal_BAF <- normal_VCF[, c("CHROM", "POS", "Sample")]
tumor_BAF  <-  tumor_VCF[, c("CHROM", "POS", "Sample")]

format_parser <- function(x){
  AD <- unlist(strsplit(x, ":"))[2]
  AD <- unlist(strsplit(AD,","))
  return(as.numeric(AD))
}

normal_AD <- lapply(normal_BAF$Sample, format_parser)
normal_AD <- setNames(do.call(rbind.data.frame, normal_AD), c("ref","baf"))

tumor_AD <- lapply(tumor_BAF$Sample, format_parser)
tumor_AD <- setNames(do.call(rbind.data.frame, tumor_AD), c("ref","baf"))

normal_BAF <- cbind(normal_BAF, normal_AD)
tumor_BAF <- cbind(tumor_BAF, tumor_AD)

# Process BAF dfs --------------------------------------

## add coverage info
normal_BAF$coverage <- normal_BAF$ref + normal_BAF$baf
tumor_BAF$coverage <- tumor_BAF$ref + tumor_BAF$baf

## keep only sites annotated in both
n_ids <- paste0(normal_BAF$CHROM, normal_BAF$POS)
t_ids <- paste0(tumor_BAF$CHROM, tumor_BAF$POS)
normal_BAF <- normal_BAF[n_ids %in% intersect(n_ids, t_ids), ]
tumor_BAF <- tumor_BAF[t_ids %in% intersect(n_ids, t_ids), ]

# Output for THeta --------------------------------------------------------
normal_Theta <- normal_BAF[, c("CHROM", "POS", "ref", "baf")]
tumor_Theta <- tumor_BAF[, c("CHROM", "POS", "ref", "baf")]

colnames(normal_Theta) <- c("#Chrm", "Pos", "Ref_Allele", "Mut_Allele")
colnames(tumor_Theta) <- c("#Chrm", "Pos", "Ref_Allele", "Mut_Allele")

normal_Theta$`#Chrm` <- sub("chr", "", normal_Theta$`#Chrm`)
tumor_Theta$`#Chrm` <- sub("chr", "", tumor_Theta$`#Chrm`)

normal_Theta$`#Chrm`[normal_Theta$`#Chrm`=="X"] <- 23
normal_Theta$`#Chrm`[normal_Theta$`#Chrm`=="Y"] <- 24

tumor_Theta$`#Chrm`[tumor_Theta$`#Chrm`=="X"] <- 23
tumor_Theta$`#Chrm`[tumor_Theta$`#Chrm`=="Y"] <- 24

dir.create("../../THetA/")
write.table(normal_Theta, file = "../../THetA/normal_SNP.txt", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(tumor_Theta, file = "../../THetA/tumor_SNP.txt", sep = "\t", row.names = F, quote = F, col.names = T)

# Output for ExomeCNV -----------------------------------------------------
normal_BAF <- normal_BAF[, c("CHROM", "POS", "coverage", "baf")]
tumor_BAF <- tumor_BAF[, c("CHROM", "POS", "coverage", "baf")]

colnames(normal_BAF) <- c("chr", "position", "coverage", "baf")
colnames(tumor_BAF) <- c("chr", "position", "coverage", "baf")

save(normal_BAF, tumor_BAF, file = "BAF_data.Rdata")
print("baf files created!")
