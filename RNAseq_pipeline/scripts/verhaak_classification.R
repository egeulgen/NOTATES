#!/usr/bin/env Rscript

#' Prepare a log2-transformed expression matrix from RSEM gene results
#'
#' Reads an RSEM `{sample}_rsem.genes.results` file and constructs a
#' single-sample expression matrix suitable for downstream analyses such as
#' ssGSEA or gene set enrichment. The function extracts TPM values, removes
#' Ensembl version numbers from gene IDs, and returns a log2(TPM + 1)
#' transformed matrix.
#'
#' @param rsem_file Character. Path to an RSEM gene results file
#'   (e.g. `{sample}_rsem.genes.results`).
#'
#' @return A numeric matrix with genes as rows and a single column
#'   (`tumour_sample`). Values are log2(TPM + 1) transformed expression.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reads the RSEM results file using `readr::read_tsv()`.
#'   \item Extracts `gene_id` and `TPM` columns.
#'   \item Removes Ensembl version numbers (e.g. `ENSG000001.5` → `ENSG000001`).
#'   \item Creates a single-column expression matrix.
#'   \item Applies log2(TPM + 1) transformation.
#' }
prepare_expression_matrix <- function(rsem_file) {
  #######################################
  # Load RSEM gene expression
  #######################################

  cat("Reading RSEM results...\n")
  rsem <- read.delim(rsem_file)
  expr <- rsem[, c("gene_id", "TPM")]

  #######################################
  # Prepare expression matrix
  #######################################

  cat("Preparing expression matrix...\n")

  expr_mat <- matrix(expr$TPM, ncol = 1)
  rownames(expr_mat) <- expr$gene_id
  colnames(expr_mat) <- "tumour_sample"

  expr_mat <- log2(expr_mat + 1)

  return(expr_mat)
}

library(msigdbr)
library(GSVA)
library(GSEABase)
library(dplyr)
library(tibble)
library(ggplot2)

#######################################
# Logging setup
#######################################

log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("Starting Verhaak classification\n")

#######################################
# Input files
#######################################

rsem_file <- snakemake@input[["rsem_file"]]


#######################################
# Output files
#######################################

out_scores_file <- snakemake@output[["scores_tsv"]]
out_scores_plot <- snakemake@output[["scores_plot"]]


#######################################
# Prepare expression matrix
#######################################
expr_mat <- prepare_expression_matrix(rsem_file)

#######################################
# Load Verhaak gene sets
#######################################

cat("Loading gene sets...\n")

cgp_gene_sets <- msigdbr(
  species = "human",
  collection = "C2",
  subcollection = "CGP"
)
verhaak_gene_sets_df <- cgp_gene_sets %>%
  filter(grepl("VERHAAK_GLIOBLASTOMA", gs_name))

gene_sets <- verhaak_gene_sets_df %>%
  mutate(
    subtype = gs_name |>
      sub("^VERHAAK_GLIOBLASTOMA_", "", x = _) |>
      sub("\\..*$", "", x = _) |>
      tolower()
  ) %>%
  group_by(subtype) %>%
  summarise(genes = list(ensembl_gene), .groups = "drop") %>%
  tibble::deframe()


#######################################
# Run ssGSEA
#######################################

cat("Running ssGSEA scoring...\n")

param <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gene_sets,
)

scores <- gsva(param)

#######################################
# Determine subtype
#######################################

scores_vec <- scores[, 1]

subtype <- names(scores_vec)[which.max(scores_vec)]

#######################################
# Output
#######################################

cat("\nVerhaak subtype scores:\n")
print(scores_vec)

cat("\nPredicted subtype:\n")
print(subtype)

#######################################
# Save results
#######################################

write.table(
  scores_vec,
  file = out_scores_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

#######################################
# Plot
#######################################

plot_df <- tibble(
  subtype = names(scores_vec),
  score = as.numeric(scores_vec)
)

g <- ggplot(plot_df, aes(x = reorder(subtype, score), y = score, fill = score))
g <- g + geom_col(width = 0.6, show.legend = FALSE)
g <- g + coord_flip()
g <- g + scale_fill_gradient(low = "lightblue", high = "steelblue")
g <- g +
  labs(
    title = "Verhaak GBM Subtype Scores",
    x = "Subtype",
    y = "ssGSEA Score"
  )
g <- g + theme_minimal(base_size = 14)

ggsave(out_scores_plot, plot = g, width = 6, height = 4)


#######################################
# Logging close
#######################################
sink(type = "message")
sink()
close(log_con)
