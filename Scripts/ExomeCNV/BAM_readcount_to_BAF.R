############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                    ExomeCNV Analysis                     #
#                BAM read-count to baf files               #
#                   Ege Ulgen, Oct 2016                    #
############################################################

# arguments from bash -----------------------------------------------------
arguments <- commandArgs(trailingOnly = T)
currentdir <- getwd()
min_cov <- arguments[1]
  
setwd(paste0(currentdir,"/ExomeCNV/baf"))

# read BAM read-count files -----------------------------------------------
## no header, find the largest row for number of columns
tmp <- readLines("./normal_bamreadcount.txt")
tmp <- strsplit(tmp, "\t")
coln <- max(sapply(tmp, length))
if(coln==10){
  cnames <-  c("chr", "position", "ref", "coverage", "eq", "A", "C", "G", "T", "N")
}else{
  cnames <-  c("chr", "position", "ref", "coverage", "eq", "A", "C", "G", "T", "N", paste0("InDel", 1:(coln - 10)))
}

normal.read_count <- read.delim("./normal_bamreadcount.txt", stringsAsFactors = F, header = F, 
                                fill = TRUE, col.names = cnames)

tmp <- readLines("./tumor_bamreadcount.txt")
tmp <- strsplit(tmp, "\t")
coln <- max(sapply(tmp, length))
if(coln==10){
  cnames <-  c("chr", "position", "ref", "coverage", "eq", "A", "C", "G", "T", "N")
}else{
  cnames <-  c("chr", "position", "ref", "coverage", "eq", "A", "C", "G", "T", "N", paste0("InDel", 1:(coln - 10)))
}

tumor.read_count <- read.delim("./tumor_bamreadcount.txt", stringsAsFactors = F, header = F, 
                                fill = TRUE, col.names = cnames)
rm(list = setdiff(ls(), c("normal.read_count", "tumor.read_count", "min_cov")))

# Process and prepare baf dataframes --------------------------------------
## minimum coverage threshold (default = 20)
normal.read_count <- normal.read_count[normal.read_count$coverage >= min_cov,]
tumor.read_count <- tumor.read_count[tumor.read_count$coverage >= min_cov,]

## Calculate B-Allele Counts
parser <- function(x){
  idx <- which(colnames(normal.read_count) == toupper(as.character(x[3])))
  tmp <- as.character(x[idx])
  ref <- as.numeric(unlist(strsplit(tmp, ":"))[2])
  baf <- as.numeric(x[4]) - ref
  return(as.numeric(baf))
}
normal.read_count$baf <- apply(normal.read_count, 1, parser)
tumor.read_count$baf <- apply(tumor.read_count, 1, parser)

## keep sites that are heterozygous in normal
baf <- normal.read_count$baf/normal.read_count$coverage
sdev <- sd(baf)
normal.read_count <- normal.read_count[baf <= 0.5 + 2*sdev & baf >= 0.5 - 2*sdev, ]

#### Create BAf dataframes and save
normal_BAF <- normal.read_count[, c("chr", "position", "coverage", "baf")]
tumor_BAF <- tumor.read_count[, c("chr", "position", "coverage", "baf")]

## keep only sites annotated in both
n_ids <- paste0(normal_BAF$chr, normal_BAF$position)
t_ids <- paste0(tumor_BAF$chr, tumor_BAF$position)
normal_BAF <- normal_BAF[n_ids %in% intersect(n_ids, t_ids),]
tumor_BAF <- tumor_BAF[t_ids %in% intersect(n_ids, t_ids),]

save(normal_BAF, tumor_BAF, file = "BAF_data.Rdata")
print("baf files created!")

# create files for THetA --------------------------------------------------
normal_SNP <- normal_BAF
tumor_SNP <- tumor_BAF

colnames(normal_SNP) <- c("chromosome", "position", "ref_allele", "mut_allele")
colnames(tumor_SNP) <- c("chromosome", "position", "ref_allele", "mut_allele")

normal_SNP$ref_allele <- normal_SNP$ref_allele - normal_SNP$mut_allele
tumor_SNP$ref_allele <- tumor_SNP$ref_allele - tumor_SNP$mut_allele

dir.create("../../THetA/")
write.table(normal_SNP, file = "../../THetA/normal_SNP.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(tumor_SNP, file = "../../THetA/tumor_SNP.txt", sep = "\t", row.names = F, quote = F, col.names = F)
