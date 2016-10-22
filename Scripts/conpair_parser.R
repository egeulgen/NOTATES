arguments <- commandArgs(trailingOnly = T)
normal_name <- arguments[1]
tumor_name <- arguments[2]

conpair <- read.delim("conpair_contam.txt", stringsAsFactors = F, header = F)

normal_cont <- unlist(strsplit(conpair[1,], ": "))[2]
tumor_cont <- unlist(strsplit(conpair[2,], ": "))[2]

normal_cont <- sub("%", "", normal_cont)
normal_cont <- as.numeric(normal_cont)/100
normal_cont <- round(normal_cont, 2)

tumor_cont <- sub("%", "", tumor_cont)
tumor_cont <- as.numeric(tumor_cont)/100
tumor_cont <- round(tumor_cont, 2)

output <- data.frame(Sample = c(normal_name, tumor_name), Contamination = c(normal_cont, tumor_cont), stringsAsFactors = F)
write.table(output, "contamination.txt", sep = "\t", quote = F, col.names = F, row.names = F)