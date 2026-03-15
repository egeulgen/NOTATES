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
