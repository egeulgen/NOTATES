##################################################
## Project: NOTATES
## Script purpose: Parse VCF to create BAF files
## Date: Oct 27, 2019
## Author: Ege Ulgen
##################################################

setwd("./ExomeCNV/baf")

# read VCF files -----------------------------------------------
cnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample")
normal_VCF <- read.delim("./normal_HQ_SNPs.vcf", header = F, col.names = cnames, comment.char = "#", stringsAsFactors = F)
tumor_VCF <- read.delim("./tumor_HQ_SNPs.vcf",  header = F, col.names = cnames, comment.char = "#", stringsAsFactors = F)

# parse and add coverages -------------------------------------------------
normal_BAF <- normal_VCF[, c("CHROM", "POS", "REF", "ALT", "Sample")]
tumor_BAF  <-  tumor_VCF[, c("CHROM", "POS", "REF", "ALT", "Sample")]

AD_format_parser <- function(x){
  AD <- unlist(strsplit(x, ":"))[2]
  AD <- unlist(strsplit(AD,","))
  return(as.numeric(AD))
}

normal_AD <- lapply(normal_BAF$Sample, AD_format_parser)
normal_AD <- setNames(do.call(rbind.data.frame, normal_AD), c("ref","baf"))

tumor_AD <- lapply(tumor_BAF$Sample, AD_format_parser)
tumor_AD <- setNames(do.call(rbind.data.frame, tumor_AD), c("ref","baf"))

normal_BAF <- cbind(normal_BAF, normal_AD)
tumor_BAF <- cbind(tumor_BAF, tumor_AD)

# Process BAF dfs ---------------------------------------------------------
## add coverage info
normal_BAF$coverage <- normal_BAF$ref + normal_BAF$baf
tumor_BAF$coverage <- tumor_BAF$ref + tumor_BAF$baf

## keep only sites annotated in both
merged_df <- merge(normal_BAF, tumor_BAF, 
                   by = c("CHROM", "POS", "REF", "ALT"),
                   suffixes = c("_N", "_T"))

normal_BAF <- merged_df[, c("CHROM", "POS", "ref_N", "baf_N", "coverage_N")]
tumor_BAF <- merged_df[c("CHROM", "POS", "ref_T", "baf_T", "coverage_T")]

# Output for THeta --------------------------------------------------------
normal_Theta <- normal_BAF[, c("CHROM", "POS", "ref_N", "baf_N")]
colnames(normal_Theta) <- c("#Chrm", "Pos", "Ref_Allele", "Mut_Allele")
tumor_Theta <- tumor_BAF[, c("CHROM", "POS", "ref_T", "baf_T")]
colnames(tumor_Theta) <- c("#Chrm", "Pos", "Ref_Allele", "Mut_Allele")

normal_Theta$`#Chrm` <- sub("chr", "", normal_Theta$`#Chrm`)
tumor_Theta$`#Chrm` <- sub("chr", "", tumor_Theta$`#Chrm`)

normal_Theta$`#Chrm`[normal_Theta$`#Chrm`=="X"] <- 23
normal_Theta$`#Chrm`[normal_Theta$`#Chrm`=="Y"] <- 24

tumor_Theta$`#Chrm`[tumor_Theta$`#Chrm`=="X"] <- 23
tumor_Theta$`#Chrm`[tumor_Theta$`#Chrm`=="Y"] <- 24

dir.create("../../THetA/")
write.table(normal_Theta, file = "../../THetA/normal_SNP.txt", 
			sep = "\t", row.names = FALSE, col.names = TRUE,
			quote = FALSE)
write.table(tumor_Theta, file = "../../THetA/tumor_SNP.txt",
			sep = "\t", row.names = FALSE, col.names = TRUE,
			quote = FALSE)

cat("\nBAF files saved for THetA!\n")

# Output for ExomeCNV -----------------------------------------------------
normal_BAF <- normal_BAF[, c("CHROM", "POS", "coverage_N", "baf_N")]
colnames(normal_BAF) <- c("chr", "position", "coverage", "baf")

tumor_BAF <- tumor_BAF[, c("CHROM", "POS", "coverage_T", "baf_T")]
colnames(tumor_BAF) <- c("chr", "position", "coverage", "baf")

save(normal_BAF, tumor_BAF, file = "BAF_data.Rdata")
cat("\nBAF files created!\n")
