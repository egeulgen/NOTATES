#!/usr/bin/env Rscript

library(quantiseqr)
library(biomaRt)
library(ggplot2)

#######################################
# Logging setup
#######################################

log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("Starting Circos Plot generation\n")


#######################################
# Inputs
#######################################

rsem_file <- snakemake@input[["rsem_file"]]
sample_id <- snakemake@params[["sample_id"]]

#######################################
# Outputs
#######################################

quanti_plot_out <- snakemake@output[["quanti_plot"]]


#######################################
# Prepare input matrix
#######################################

cat("Preparing input matrix...\n")

rsem_result_df <- read.table(rsem_file, header = TRUE, as.is = TRUE)
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rsem_result_df$gene_id,
  mart = ensembl
)
rsem_result_df$ensembl_gene_id <- rsem_result_df$gene_id
tpm_mapped <- merge(
  mapping,
  rsem_result_df[, c("ensembl_gene_id", "TPM")],
  by = "ensembl_gene_id"
)
tpm_mapped <- tpm_mapped[tpm_mapped$hgnc_symbol != "", ]
tpm_mapped <- tpm_mapped[!duplicated(tpm_mapped$hgnc_symbol), ]
rownames(tpm_mapped) <- tpm_mapped$hgnc_symbol
tpm_mapped <- tpm_mapped[, "TPM", drop = FALSE]


cat("Running quanTIseq...\n")

result <- quantiseqr::run_quantiseq(
  expression_data = tpm_mapped,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

result$Sample <- sample_id
rownames(result) <- sample_id
g <- quantiplot(result)
ggsave(quanti_plot_out, g, width = 15, height = 5)
