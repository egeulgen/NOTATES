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

cat("Starting enrichment analyses\n")

#######################################
# Input files
#######################################

rsem_file <- snakemake@input[["rsem_file"]]

#######################################
# Output files
#######################################
collections <- snakemake@params[["collections"]]

# Named list of output TSV files per collection
out_scores_files <- list()
out_scores_plots <- list()
for (col in collections) {
  out_scores_files[[col]] <- snakemake@output[[paste0(col, "_scores_tsv")]]
  out_scores_plots[[col]] <- snakemake@output[[paste0(col, "_scores_pdf")]]
}

#######################################
# Prepare expression matrix
#######################################
expr_mat <- prepare_expression_matrix(rsem_file)

#######################################
# Define collections
#######################################
for (col in collections) {
  cat("\nProcessing collection:", col, "\n")

  # Load gene sets for this collection
  if (col == "KEGG") {
    gene_sets_df <- msigdbr(
      species = "human",
      collection = "C2",
      subcollection = "CP",
      subcategory = "KEGG_LEGACY"
    )
  } else {
    gene_sets_df <- msigdbr(species = "human", collection = col)
  }

  gene_sets <- gene_sets_df %>%
    group_by(gs_name) %>%
    summarise(genes = list(ensembl_gene), .groups = "drop") %>%
    tibble::deframe()

  # Run ssGSEA
  param <- ssgseaParam(exprData = expr_mat, geneSets = gene_sets)
  scores <- gsva(param)

  # Convert to named vector per sample (only first column if single sample)
  scores_vec <- scores[, 1]

  # Save all scores sorted decreasingly
  out_file <- out_scores_files[[col]]
  scores_sorted <- sort(scores_vec, decreasing = TRUE)
  write.table(
    scores_sorted,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  # Top 25 entries
  top25 <- sort(scores_vec, decreasing = TRUE)[1:25]
  cat("\nTop 25 enriched gene sets for collection", col, ":\n")
  print(top25)

  # Plot top 25 using ggplot2
  top10_df <- tibble(
    gs_name = names(top25),
    score = as.numeric(top25)
  )

  g <- ggplot(
    top10_df,
    aes(x = reorder(gs_name, score), y = score, fill = score)
  )
  g <- g + geom_col(width = 0.6, show.legend = FALSE)
  g <- g + coord_flip()
  g <- g + scale_fill_gradient(low = "lightblue", high = "steelblue")
  g <- g +
    labs(
      title = paste("Top 25 ssGSEA scores - Collection", col),
      x = "Gene set",
      y = "ssGSEA score"
    )
  g <- g + theme_minimal(base_size = 12)

  ggsave(out_scores_plots[[col]], plot = g, width = 16, height = 10)
}

cat("\nAll collections processed.\n")

#######################################
# Logging close
#######################################
sink(type = "message")
sink()
close(log_con)
