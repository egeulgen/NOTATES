##################################################
## Project: NOTATES
## Script purpose: Annotation of genes within segments 
## and the cytobands containing the segments
## Date: Sep 22, 2020
## Author: Ege Ulgen
##################################################

# Installations etc. ------------------------------------------------------
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GenomicRanges))

# SCNA Annotation ---------------------------------------------------------

#' Annotate Genes within Segments (hg38)
#'
#' @param segment_df data frame of segments (first 3 columns must be segment's chr, start and end)
#'
#' @return a list containing 2 objects: \describe{
#'   \item{Gene_segments}{list of segments per gene, each element is a vector of segment intervals}
#'   \item{Segment_genes}{list of genes per segment, each element containg the genes in the given segment}
#' }
annotate_genes <- function(segment_df) {
  segment_gr <- makeGRangesFromDataFrame(segment_df)
  hg38_genes_gr <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  # find regions overlapping genes
  overlapGenes <- findOverlaps(query = hg38_genes_gr, 
                               subject = segment_gr, 
                               type = "within",
                               ignore.strand = TRUE)
  
  #####  list genes within each segment
  tmp <- factor(subjectHits(overlapGenes),
                levels = seq_len(subjectLength(overlapGenes)))
  res <- splitAsList(mcols(hg38_genes_gr)[["gene_id"]][queryHits(overlapGenes)], tmp)
  res <- as.list(res)
  
  all_syms_tbl <- as.data.frame(org.Hs.egSYMBOL)
  segment_genes <- lapply(res, function(x) {
    if (length(x) == 0)
      return(NULL)
    # convert to gene symbol
    syms <- all_syms_tbl$symbol[match(x, all_syms_tbl$gene_id)]
    syms <- syms[!is.na(syms)]
    return(syms)
  })
  
  ##### list overlapping segments for each gene
  segment_gr$seg_idx <- seq_len(nrow(segment_df))
  
  tmp <- factor(queryHits(overlapGenes),
                levels = seq_len(queryLength(overlapGenes)))
  res <- splitAsList(mcols(segment_gr)[["seg_idx"]][subjectHits(overlapGenes)], tmp)
  res <- as.list(res)
  
  # exclude genes with no segments
  cnts <- vapply(res, length, 1L)
  res <- res[cnts != 0]
  
  # fetch gene ids using idx
  names(res) <- mcols(hg38_genes_gr)[["gene_id"]][as.numeric(names(res))]
  
  # convert ids to symbols
  tmp <- all_syms_tbl$symbol[match(names(res), all_syms_tbl$gene_id)]
  res <- res[!is.na(tmp)]
  names(res) <- all_syms_tbl$symbol[match(names(res), all_syms_tbl$gene_id)]
  
  gene_segments <- lapply(res, function(x) {
    tmp <- segment_df[x, ]
    seg_pos_vec <- apply(tmp, 1, function(y) paste0(y[1], ":", y[2], "-", y[3]))
    seg_pos_vec <- unname(seg_pos_vec)
    return(seg_pos_vec)
  })
  
  return(list(Gene_segments = gene_segments, 
              Segment_genes = segment_genes))
}


#' Annotate Segments Overlapping Cytobands (hg38)
#'
#' @param segment_df data frame of segments (first 3 columns must be segment's chr, start and end)
#' @param cytobands_df data frame of cytobands (can be obtained through UCSC)
#' 
#' @return a list containing 2 objects: \describe{
#'   \item{Cytb_segments}{list of segments per cytoband, each element is a vector of segment intervals}
#'   \item{Segment_cytbs}{list of cytobands overlapped by the given segment, each element is a vector of cytobands for the segment}
#' }
annotate_cytb <- function(segment_df, cytobands_df) {
  
  ## Segments GRanges object
  segment_gr <- makeGRangesFromDataFrame(segment_df)
  segment_gr$seg_idx <- seq_len(nrow(segment_df))
  
  ## Cytobands GRanges object
  cytb_gr <- makeGRangesFromDataFrame(cytobands_df,
                                      starts.in.df.are.0based = TRUE,
                                      keep.extra.columns = TRUE)
  
  # find regions overlapping genes
  olap <- findOverlaps(query = cytb_gr, 
                       subject = segment_gr, 
                       type = "within",
                       ignore.strand = TRUE)
  
  ##### list cytoband(s) within each segment
  tmp <- factor(subjectHits(olap),
                levels = seq_len(subjectLength(olap)))
  res <- splitAsList(mcols(cytb_gr)[["Cytb_name"]][queryHits(olap)], tmp)
  res <- as.list(res)
  
  segment_cytbs <- lapply(res, function(x) {
    if (length(x) == 0)
      return(NULL)
    return(x)
  })
  
  ##### list overlapping segments for each cytoband
  tmp <- factor(queryHits(olap),
                levels = seq_len(queryLength(olap)))
  res <- splitAsList(mcols(segment_gr)[["seg_idx"]][subjectHits(olap)], tmp)
  res <- as.list(res)
  
  # exclude genes with no segments
  cnts <- vapply(res, length, 1L)
  res <- res[cnts != 0]
  
  # fetch gene ids using idx
  names(res) <- mcols(cytb_gr)[["Cytb_name"]][as.numeric(names(res))]
  
  cytb_segments <- lapply(res, function(x) {
    tmp <- segment_df[x, ]
    seg_pos_vec <- apply(tmp, 1, function(y) paste0(y[1], ":", y[2], "-", y[3]))
    seg_pos_vec <- unname(seg_pos_vec)
    return(seg_pos_vec)
  })
  
  return(list(Cytb_segments = cytb_segments, 
              Segment_cytbs = segment_cytbs))
}

# LOH annotation ----------------------------------------------------------

#' LOH - Annotate Genes overlapping Segments (hg38)
#'
#' @param segment_df data frame of segments (first 3 columns must be segment's chr, start and end)
#'
#' @return a list containing 2 objects: \describe{
#'   \item{Gene_segments}{list of segments per gene, each element is a vector of segment intervals}
#'   \item{Segment_genes}{list of genes per segment, each element containg the genes in the given segment}
#' }
loh_annotate_genes <- function(segment_df) {
  segment_gr <- makeGRangesFromDataFrame(segment_df)
  hg38_genes_gr <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  # find regions overlapping genes
  overlapGenes <- findOverlaps(query = hg38_genes_gr, 
                               subject = segment_gr, 
                               type = "any",
                               ignore.strand = TRUE)
  
  #####  list genes within each segment
  tmp <- factor(subjectHits(overlapGenes),
                levels = seq_len(subjectLength(overlapGenes)))
  res <- splitAsList(mcols(hg38_genes_gr)[["gene_id"]][queryHits(overlapGenes)], tmp)
  res <- as.list(res)
  
  all_syms_tbl <- as.data.frame(org.Hs.egSYMBOL)
  segment_genes <- lapply(res, function(x) {
    if (length(x) == 0)
      return(NULL)
    # convert to gene symbol
    syms <- all_syms_tbl$symbol[match(x, all_syms_tbl$gene_id)]
    syms <- syms[!is.na(syms)]
    return(syms)
  })
  
  ##### list overlapping segments for each gene
  segment_gr$seg_idx <- seq_len(nrow(segment_df))
  
  tmp <- factor(queryHits(overlapGenes),
                levels = seq_len(queryLength(overlapGenes)))
  res <- splitAsList(mcols(segment_gr)[["seg_idx"]][subjectHits(overlapGenes)], tmp)
  res <- as.list(res)
  
  # exclude genes with no segments
  cnts <- vapply(res, length, 1L)
  res <- res[cnts != 0]
  
  # fetch gene ids using idx
  names(res) <- mcols(hg38_genes_gr)[["gene_id"]][as.numeric(names(res))]
  
  # convert ids to symbols
  tmp <- all_syms_tbl$symbol[match(names(res), all_syms_tbl$gene_id)]
  res <- res[!is.na(tmp)]
  names(res) <- all_syms_tbl$symbol[match(names(res), all_syms_tbl$gene_id)]
  
  gene_segments <- lapply(res, function(x) {
    tmp <- segment_df[x, ]
    seg_pos_vec <- apply(tmp, 1, function(y) paste0(y[1], ":", y[2], "-", y[3]))
    seg_pos_vec <- unname(seg_pos_vec)
    return(seg_pos_vec)
  })
  
  return(list(Gene_segments = gene_segments, 
              Segment_genes = segment_genes))
}


#' LOH - Annotate Segments Overlapping Cytobands (hg38)
#'
#' @param segment_df data frame of segments (first 3 columns must be segment's chr, start and end)
#' @param cytobands_df data frame of cytobands (can be obtained through UCSC)
#' 
#' @return a list containing 2 objects: \describe{
#'   \item{Cytb_segments}{list of segments per cytoband, each element is a vector of segment intervals}
#'   \item{Segment_cytbs}{list of cytobands overlapped by the given segment, each element is a vector of cytobands for the segment}
#' }
loh_annotate_cytb <- function(segment_df, cytobands_df) {
  
  ## Segments GRanges object
  segment_gr <- makeGRangesFromDataFrame(segment_df)
  segment_gr$seg_idx <- seq_len(nrow(segment_df))
  
  ## Cytobands GRanges object
  cytb_gr <- makeGRangesFromDataFrame(cytobands_df,
                                      starts.in.df.are.0based = TRUE,
                                      keep.extra.columns = TRUE)
  
  # find regions overlapping genes
  olap <- findOverlaps(query = cytb_gr, 
                       subject = segment_gr, 
                       type = "any",
                       ignore.strand = TRUE)
  
  ##### list cytoband(s) overlapping each segment
  tmp <- factor(subjectHits(olap),
                levels = seq_len(subjectLength(olap)))
  res <- splitAsList(mcols(cytb_gr)[["Cytb_name"]][queryHits(olap)], tmp)
  res <- as.list(res)
  
  segment_cytbs <- lapply(res, function(x) {
    if (length(x) == 0)
      return(NULL)
    return(x)
  })
  
  ##### list overlapping segments for each cytoband
  tmp <- factor(queryHits(olap),
                levels = seq_len(queryLength(olap)))
  res <- splitAsList(mcols(segment_gr)[["seg_idx"]][subjectHits(olap)], tmp)
  res <- as.list(res)
  
  # exclude genes with no segments
  cnts <- vapply(res, length, 1L)
  res <- res[cnts != 0]
  
  # fetch gene ids using idx
  names(res) <- mcols(cytb_gr)[["Cytb_name"]][as.numeric(names(res))]
  
  cytb_segments <- lapply(res, function(x) {
    tmp <- segment_df[x, ]
    seg_pos_vec <- apply(tmp, 1, function(y) paste0(y[1], ":", y[2], "-", y[3]))
    seg_pos_vec <- unname(seg_pos_vec)
    return(seg_pos_vec)
  })
  
  return(list(Cytb_segments = cytb_segments, 
              Segment_cytbs = segment_cytbs))
}


# Common ------------------------------------------------------------------

#' Shorten Cytoband Annotation
#'
#' @param cyt vector of cytoband annotations (character)
#' @param cytobands_df data frame of cytobands (can be obtained through UCSC)
#'
#' @return string containing shortened cytoband annotation. e.g.
#' c("chr1p36.33", "chr1p36.32", "chr1p36.31", "chr1p36.23") >> "chr1p36.33-chr1p36.23"
shorten_cyt_anno <- function(cyt, cytobands_df)
{
  if(length(cyt) == 1) ## single cytoband
    return(cyt)
  if(is.null(cyt))
    return(NA)
  
  tmp <- subset(cytobands_df, Chr == sub("[p-q].*", "", cyt[1]))
  
  tmp_P <- tmp$Cytb_name[grepl("p", tmp$Cytb_name)]
  tmp_P_start <- tmp_P[1]
  tmp_P_end <- tmp_P[length(tmp_P)]
  
  tmp_Q <- tmp$Cytb_name[grepl("q", tmp$Cytb_name)]
  tmp_Q_start <- tmp_Q[1]
  tmp_Q_end <- tmp_Q[length(tmp_Q)]
  
  ## whole chromosome
  if (cyt[1] == tmp_P_start & cyt[length(cyt)] == tmp_Q_end) 
    return(sub("[p-q].*", "", cyt[1]))
  ## whole short arm
  if (cyt[1] == tmp_P_start & cyt[length(cyt)] == tmp_P_end) 
    return(paste0(sub("[p-q].*", "", cyt[1]), "p"))
  ## whole long arm
  if (cyt[1] == tmp_Q_start & cyt[length(cyt)] == tmp_Q_end) 
    return(paste0(sub("[p-q].*", "", cyt[1]), "q"))
  ## multiple segments
  return(paste0(cyt[1], "-", cyt[length(cyt)]))
}
