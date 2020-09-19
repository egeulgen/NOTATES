##################################################
## Project: NOTATES
## Script purpose: Prepare WES QC metrics table
## Date: Oct 27, 2019
## Author: Ege Ulgen
##################################################

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
normal_name <- args[1]
tumor_name <- args[2]

# Create QC tables --------------------------------------------------------
create_QC_table <- function(sample_name, type, r_type = "Paired-End") {
  #### Number of lanes
  numlanes <- readLines(file.path(sample_name, "lanes.txt"))
  numlanes <- length(numlanes)
  
  #### Alignment Metrics
  alignment_metrics <- read.delim(file.path(sample_name, "QC/alignment_summary_metrics.txt"),
                                  comment.char = "#", nrows = 3)
  algn_summ_keep <- c("MEAN_READ_LENGTH", # Mean read length
                      "TOTAL_READS", # The total number of reads
                      "PCT_PF_READS", # The percentage of reads that are PF where PF is defined as passing Illumina's filter
                      "PCT_PF_READS_ALIGNED", # The percentage of PF reads that aligned to the reference sequence
                      "PCT_READS_ALIGNED_IN_PAIRS") # The percentage of reads whose mate pair was also aligned to the reference
                      
  alignment_metrics <- alignment_metrics[alignment_metrics$CATEGORY == "PAIR", 
                                         algn_summ_keep]

  ### Coverage metrics
  fname <- paste0(type, ".coverage.sample_summary")
  coverage_metrics <- read.delim(file.path("./ExomeCNV/DepthOfCoverage", fname))
  
  cov_keep <- c("mean", "X._bases_above_1", "X._bases_above_5", "X._bases_above_10", 
                "X._bases_above_25", "X._bases_above_50", "X._bases_above_100")
  coverage_metrics <- coverage_metrics[coverage_metrics$sample_id == sample_name, cov_keep]
  
  # combine
  sample_QC_tbl <- data.frame(Number_of_Lanes = numlanes,
                              Read_type = r_type,
                              alignment_metrics,
                              coverage_metrics)
  
  # format
  sample_QC_tbl$TOTAL_READS <- round(sample_QC_tbl$TOTAL_READS / 1e6, 2)
  sample_QC_tbl$PCT_PF_READS <- round(sample_QC_tbl$PCT_PF_READS * 100, 2)
  sample_QC_tbl$PCT_PF_READS_ALIGNED <- round(sample_QC_tbl$PCT_PF_READS_ALIGNED * 100, 2)
  sample_QC_tbl$PCT_READS_ALIGNED_IN_PAIRS <- round(sample_QC_tbl$PCT_READS_ALIGNED_IN_PAIRS * 100, 2)
  
  return(sample_QC_tbl)
}

normal_QC_table <- create_QC_table(normal_name, "normal")
tumor_QC_table <- create_QC_table(tumor_name, "tumor")
  
# Merge and format --------------------------------------------------------
QC_table <- rbind(normal_QC_table, tumor_QC_table)
# format
QC_table$PCT_PF_READS <- paste0(QC_table$PCT_PF_READS, "%")
QC_table$PCT_PF_READS_ALIGNED <- paste0(QC_table$PCT_PF_READS_ALIGNED, "%")
QC_table$PCT_READS_ALIGNED_IN_PAIRS <- paste0(QC_table$PCT_READS_ALIGNED_IN_PAIRS, "%")

QC_table$X._bases_above_1 <- paste0(QC_table$X._bases_above_1, "%")
QC_table$X._bases_above_5 <- paste0(QC_table$X._bases_above_5, "%")
QC_table$X._bases_above_10 <- paste0(QC_table$X._bases_above_10, "%")
QC_table$X._bases_above_25 <- paste0(QC_table$X._bases_above_25, "%")
QC_table$X._bases_above_50<- paste0(QC_table$X._bases_above_50, "%")
QC_table$X._bases_above_100 <- paste0(QC_table$X._bases_above_100, "%")

# transpose
QC_table <- as.data.frame(t(QC_table))


colnames(QC_table) <- c("Normal", "Tumor")
rownames(QC_table) <- c("Number of Lanes","Read Type", "Read Length",
                              "Total Number of Reads(in millions)",
                              "PF Reads*", "Aligned PF Reads**", "PE Aligned***",
                              "Mean Coverage", "1X", "5X", "10X", "25X", "50X", "100X")
QC_table <- apply(QC_table, 2, trimws)

write.csv(QC_table, "QC_table.csv", quote = F)
