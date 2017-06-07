############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                     Quality Metrics                      #
#                   Ege Ulgen, Apr 2017                    #
############################################################

dir.create("./Quality_Metrics/")
setwd("./Quality_Metrics/")

args <- commandArgs(trailingOnly=TRUE)
normal_name <- args[1]
tumor_name <- args[2]

# Normal ------------------------------------------------------------------
numlanes <- list.files(paste0("../", normal_name, "/QC/"))
numlanes <- numlanes[grepl("^L", numlanes)]
numlanes <- unique(gsub("_.*","",numlanes))
numlanes <- length(numlanes)

#### Alignment Metrics
alignment_metrics <- read.delim(paste0("../", normal_name,"/QC/alignment_summary_metrics.txt"), 
                                stringsAsFactors = F, comment.char = "#")
## The total number of reads
total_reads <- alignment_metrics$TOTAL_READS[3]
## The percentage of reads that are PF where PF is defined as passing Illumina's filter
pct_pf_reads <- alignment_metrics$PCT_PF_READS[3]
## The percentage of PF reads that aligned to the reference sequence
pct_mapped <- alignment_metrics$PCT_PF_READS_ALIGNED[3]
## The percentage of reads whose mate pair was also aligned to the reference
pct_PE_mapped <- alignment_metrics$PCT_READS_ALIGNED_IN_PAIRS[3]
## Read length
read_length <- alignment_metrics$MEAN_READ_LENGTH[3]

### Coverage metrics
coverage_metrics <- read.delim(paste0("../", normal_name,"/QC/primary_target_coverage.sample_summary"),
                               stringsAsFactors = F, comment.char = "#")

mean_coverage <- coverage_metrics$mean[1]
bases_1X <- coverage_metrics$X._bases_above_1[1]
bases_5X <- coverage_metrics$X._bases_above_5[1]
bases_10X <- coverage_metrics$X._bases_above_10[1]
bases_25X <- coverage_metrics$X._bases_above_25[1]
bases_50X <- coverage_metrics$X._bases_above_50[1]
bases_100X <- coverage_metrics$X._bases_above_100[1]

normal_QC_table <- data.frame(Number_of_Lanes = numlanes, Read_type = "Paired-End", Read_length = read_length,
                              Total_Reads = total_reads, Pct_PF_reads = pct_pf_reads, Pct_mapped = pct_mapped,
                              Pct_PE_mapped = pct_PE_mapped, Mean_cov = mean_coverage, 
                              X1 = bases_1X, X5 = bases_5X, X10 = bases_10X, 
                              X25 = bases_25X, X50 = bases_50X, X100 = bases_100X)

normal_QC_table$Total_Reads <- normal_QC_table$Total_Reads/10^6
normal_QC_table$Pct_PF_reads <- round(normal_QC_table$Pct_PF_reads*100, 2)
normal_QC_table$Pct_mapped <- round(normal_QC_table$Pct_mapped*100, 2)
normal_QC_table$Pct_PE_mapped <- round(normal_QC_table$Pct_PE_mapped*100, 2)
normal_QC_table <- t(normal_QC_table)

# Tumor ------------------------------------------------------------------
numlanes <- list.files(paste0("../", tumor_name, "/QC/"))
numlanes <- numlanes[grepl("^L", numlanes)]
numlanes <- unique(gsub("_.*","",numlanes))
numlanes <- length(numlanes)

#### Alignment Metrics
alignment_metrics <- read.delim(paste0("../", tumor_name,"/QC/alignment_summary_metrics.txt"), 
                                stringsAsFactors = F, comment.char = "#")
## The total number of reads
total_reads <- alignment_metrics$TOTAL_READS[3]
## The percentage of reads that are PF where PF is defined as passing Illumina's filter
pct_pf_reads <- alignment_metrics$PCT_PF_READS[3]
## The percentage of PF reads that aligned to the reference sequence
pct_mapped <- alignment_metrics$PCT_PF_READS_ALIGNED[3]
## The percentage of reads whose mate pair was also aligned to the reference
pct_PE_mapped <- alignment_metrics$PCT_READS_ALIGNED_IN_PAIRS[3]
## Read length
read_length <- alignment_metrics$MEAN_READ_LENGTH[3]

### Coverage metrics
coverage_metrics <- read.delim(paste0("../", tumor_name,"/QC/primary_target_coverage.sample_summary"),
                               stringsAsFactors = F, comment.char = "#")

mean_coverage <- coverage_metrics$mean[1]
bases_1X <- coverage_metrics$X._bases_above_1[1]
bases_5X <- coverage_metrics$X._bases_above_5[1]
bases_10X <- coverage_metrics$X._bases_above_10[1]
bases_25X <- coverage_metrics$X._bases_above_25[1]
bases_50X <- coverage_metrics$X._bases_above_50[1]
bases_100X <- coverage_metrics$X._bases_above_100[1]

tumor_QC_table <- data.frame(Number_of_Lanes = numlanes, Read_type = "Paired-End", Read_length = read_length,
                              Total_Reads = total_reads, Pct_PF_reads = pct_pf_reads, Pct_mapped = pct_mapped,
                              Pct_PE_mapped = pct_PE_mapped, Mean_cov = mean_coverage, 
                              X1 = bases_1X, X5 = bases_5X, X10 = bases_10X, 
                              X25 = bases_25X, X50 = bases_50X, X100 = bases_100X)

tumor_QC_table$Total_Reads <- tumor_QC_table$Total_Reads/10^6
tumor_QC_table$Pct_PF_reads <- round(tumor_QC_table$Pct_PF_reads*100, 2)
tumor_QC_table$Pct_mapped <- round(tumor_QC_table$Pct_mapped*100, 2)
tumor_QC_table$Pct_PE_mapped <- round(tumor_QC_table$Pct_PE_mapped*100, 2)
tumor_QC_table <- t(tumor_QC_table)

# Merge -------------------------------------------------------------------
QC_table <- cbind(normal_QC_table,tumor_QC_table)
colnames(QC_table) <- c("Normal", "Tumor")
rownames(QC_table) <- c("Number of Lanes","Read Type", "Read Length",
                              "Total Number of Reads(in millions)",
                              "PF Reads* (%)", "Mapped** (%)", "PE Mapped*** (%)",
                              "Mean Coverage", "1X (%)", "5X (%)", "10X (%)", "25X (%)", "50X (%)", "100X (%)")

write.csv(QC_table, "QC_table.csv", quote = F)
