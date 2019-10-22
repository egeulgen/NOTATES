# load and install required libraries -------------------------------------
if(!require(MSIseq))
  install.packages("MSIseq")
library(MSIseq)

arguments <- commandArgs(trailingOnly = T)
patientID <- arguments[1]
exome_len <- as.numeric(arguments[2])
script_dir <- arguments[3]

## Input MAF
all_muts <- read.delim("./Oncotator/annotated.sSNVs.tsv", comment.char = "#")
all_muts <- subset(all_muts, 
                   alt_allele_seen == "True")
all_muts$Chrom <- paste0("chr", all_muts$Chromosome)
input <- all_muts[,c("Chrom", "Start_position", "End_position", "Variant_Type", "Tumor_Sample_Barcode")]
colnames(input) <- sub("position", "Position", colnames(input))
input$Tumor_Sample_Barcode <- patientID

input <- input[input$Variant_Type %in% c("SNP", "INS", "DEL"), ]

## Sequence Lengths
seq_lengths <- data.frame(Tumor_Sample_Barcode = patientID,
                          Sequence_Length = exome_len)

## hg19 repeats
# url <- "http://steverozen.net/data/Hg19repeats.rda"
# file <- paste0("~/NOTATES/Scripts/", basename(url))
# download.file(url, file)
load(paste0(script_dir, "/Hg19repeats.rda"))

dir.create("MSIseq_out")

# Input processing and classification
cat("## Computing input variables\n")
test_mutationNum <- Compute.input.variables(input, repeats = Hg19repeats, seq.len = seq_lengths)
write.csv(test_mutationNum, "MSIseq_out/Input_vars.csv")

cat("## Classifying the tumor as MSI-H or non-MSI-H\n")
result <- MSIseq.classify(mutationNum = test_mutationNum)
write.csv(result, "MSIseq_out/result.csv")

cat(paste0("The tumor is ", result$MSI_status, "\n"))
