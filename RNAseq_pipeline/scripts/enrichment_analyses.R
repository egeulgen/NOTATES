#!/usr/bin/env Rscript

source(file.path("scripts", "prepare_expression_matrix.R"))

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
collections <- c("H", "C6", "C7", "KEGG")

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

  # Top 10 entries
  top10 <- sort(scores_vec, decreasing = TRUE)[1:10]
  cat("\nTop 10 enriched gene sets for collection", col, ":\n")
  print(top10)

  # Plot top 10 using ggplot2
  top10_df <- tibble(
    gs_name = names(top10),
    score = as.numeric(top10)
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
      title = paste("Top 10 ssGSEA scores - Collection", col),
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
