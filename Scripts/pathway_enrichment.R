##################################################
## Project: NOTATES
## Script purpose: Script for performing KEGG pathway
## enrichment analysis (ORA) using high-confidence genes
## Date: Oct 27, 2019
## Author: Ege Ulgen
##################################################

# dir for data sources (same as the directory where script is located)
initial_options <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", initial_options[grep("--file=", initial_options)])
script_dir <- dirname(script_name)

# Install package(s) if necessary -----------------------------------------
if(!suppressPackageStartupMessages(require(KEGGREST))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("KEGGREST")
  suppressPackageStartupMessages(library(KEGGREST))
}

if (!suppressPackageStartupMessages(require(Homo.sapiens))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Homo.sapiens")
  suppressPackageStartupMessages(library(Homo.sapiens))
}

if(!suppressPackageStartupMessages(require(pathview))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("pathview")
  suppressPackageStartupMessages(library(pathview))
}

# Extract coding genes with mutations -------------------------------------
somatic_SNVs <- read.delim("./Oncotator/annotated.sSNVs.tsv", stringsAsFactors=F, comment.char="#")

# Subsetting for (MuTect's default) HQ filters
somatic_SNVs <- subset(somatic_SNVs, 
                       alt_allele_seen=="True" &
                       short_tandem_repeat_membership == "False")

# Keep only coding mutations
somatic_SNVs <- somatic_SNVs[!somatic_SNVs$Variant_Classification %in% c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", 
                                                                         "IGR", "Intron", "lincRNA", "RNA"),]

genes_df <- data.frame(Gene = unique(somatic_SNVs$Hugo_Symbol),
                                Value = 1, stringsAsFactors = F)

# Exract coding genes with CNV --------------------------------------------
cnv_df <- read.csv("NOTATES/SCNA/all_genes.csv", stringsAsFactors = F)
cnv_df <- cnv_df[cnv_df$ratio <= 0.5 | cnv_df$ratio >= 1.5,]

if( anyDuplicated(cnv_df$Gene) != 0 ) {
  dups <- unique(cnv_df$Gene[duplicated(cnv_df$Gene)])
  for(dup in dups) {
    cnv_df$ratio[cnv_df$Gene == dup] <- mean(cnv_df$ratio[cnv_df$Gene == dup])
  }
}
cnv_df <- cnv_df[!duplicated(cnv_df$Gene),]

cnv_df$logR <- log2(cnv_df$ratio)
cnv_df <- cnv_df[,c("Gene", "logR")]

final_df <- merge(genes_df, cnv_df, by = "Gene", all = TRUE)
final_df$Value[is.na(final_df$Value)] <- 0
final_df$logR[is.na(final_df$logR)] <- 0

genes_of_interest <- unique(final_df$Gene)

# Enrichment -----------------------------------
# created named list, eg:  path:map00010: "Glycolysis / Gluconeogenesis"
pathways_list <- keggList("pathway", "hsa")

# make them into KEGG-style human pathway identifiers
pathway_codes <- sub("path:", "", names(pathways_list))
pathways_list <- sub(" - Homo sapiens \\(human\\)", "", pathways_list)

# 
# # subsetting by c(TRUE, FALSE) -- which repeats
# # as many times as needed, sorts through some
# # unexpected packaging of geneIDs in the GENE element
# # of each pw[[n]]
# kegg_gene_sets <- sapply(pathway_codes, function(pwid){
#   pw <- keggGet(pwid)
#   pw <- pw[[1]]$GENE[c(F, T)]
#   pw <- sub(";.+", "", pw)
#   pw <- pw[grep("^[A-Za-z0-9_-]+(\\@)?$", pw)] ## removing mistaken lines
#   pw <- unique(pw)
#   return(pw)
# })
# saveRDS(kegg_gene_sets, file.path(script_dir, "kegg_gene_sets.RDS"))

kegg_gene_sets <- readRDS(file.path(script_dir, "kegg_gene_sets.RDS"))

# exclude very smmall and very large pw gene sets
cnts <- vapply(kegg_gene_sets, length, 1)
kegg_gene_sets <- kegg_gene_sets[cnts >= 10 & cnts <= 300]

# universal gene set
all_genes <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tmp <- as.data.frame(org.Hs.egSYMBOL)
all_genes <- tmp$symbol[match(all_genes$gene_id, tmp$gene_id)]
all_genes <- all_genes[!is.na(all_genes)]
all_genes <- unique(all_genes)

hyperg_test <- function(pw_genes, chosen_genes, all_genes) {
  pw_genes_selected <- length(intersect(chosen_genes, pw_genes))
  pw_genes_in_pool <- length(pw_genes)
  tot_genes_in_pool <- length(all_genes)
  non_pw_genes_in_pool <- tot_genes_in_pool - pw_genes_in_pool
  num_selected_genes <- length(chosen_genes)
  
  stats::phyper(pw_genes_selected - 1, pw_genes_in_pool,
                non_pw_genes_in_pool, num_selected_genes,
                lower.tail = FALSE)
}

enrichment_res <- sapply(kegg_gene_sets, hyperg_test, genes_of_interest, all_genes)
enrichment_res <- sort(enrichment_res)
enrichment_res <- data.frame(Pathway = pathways_list[match(paste0("path:", names(enrichment_res)), 
                                                           names(pathways_list))],
                             p = enrichment_res,
                             FDR_p = p.adjust(enrichment_res, method = "BH"), 
                             row.names = names(enrichment_res))
enrichment_res <- enrichment_res[enrichment_res$FDR_p < 0.05,, drop = F]

# Visualization -----------------------------------------------------------
rownames(final_df) <- final_df$Gene
final_df <- final_df[,-1]
colnames(final_df) <- c("Mutation", "CNV")

dir.create("pathways")
setwd("pathways/")

### Pathways in Cancer
suppressMessages(pathview(gene.data = final_df, gene.idtype = "SYMBOL",
                          pathway.id = "hsa05200", species = "hsa", out.suffix = "CANCER", 
                          keys.align = "y", kegg.native = TRUE, key.pos = "topright", same.layer = FALSE, silent = TRUE))

### Enriched Pathways 
enrichment_res$Somatic_Mutation <- ""
enrichment_res$SCNA_down <- ""
enrichment_res$SCNA_up <- ""

cnv_up <- cnv_df$Gene[cnv_df$logR > 0 ]
cnv_down <- cnv_df$Gene[cnv_df$logR < 0 ]

for(i in 1:nrow(enrichment_res))
{
  path_suffix <- gsub("\\/", "_", enrichment_res$Pathway[i])
  tmp <- suppressMessages(pathview(gene.data = final_df, gene.idtype = "SYMBOL",
                                   pathway.id = rownames(enrichment_res)[i], species = "hsa", out.suffix = path_suffix, 
                                   keys.align = "y", kegg.native = TRUE, key.pos = "topright", same.layer = FALSE, silent = TRUE))
  if (is.list(tmp)) {
    
    tmp <- tmp$plot.data.gene$all.mapped
    tmp <- tmp[tmp!=""]
    tmp <- unlist(strsplit(tmp, split = ","))
    tmp <- unique(tmp)
    tmp <- select(org.Hs.eg.db, tmp, "SYMBOL", "ENTREZID")[,2]
    
    enrichment_res$Somatic_Mutation[i] <- paste(tmp[tmp %in% genes_df$Gene], collapse = ", ")
    enrichment_res$SCNA_down[i] <- paste(tmp[tmp %in% cnv_down], collapse = ", ")
    enrichment_res$SCNA_up[i] <- paste(tmp[tmp %in% cnv_up], collapse = ", ")
  }
}

write.csv(enrichment_res, "enrichment_results.csv", row.names = T)

setwd("..")