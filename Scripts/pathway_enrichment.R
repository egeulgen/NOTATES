
# Install package(s) if necessary -----------------------------------------
if(!"KEGGREST" %in% installed.packages())
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("KEGGREST")
}

if(!"org.Hs.eg.db" %in% installed.packages())
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
}

if(!"Category" %in% installed.packages())
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("Category")
}


if(!"pathview" %in% installed.packages())
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("pathview")
}

# Extract coding genes with mutations -------------------------------------
somatic_SNVs <- read.delim("./Oncotator/annotated.sSNVs.tsv", stringsAsFactors=F, comment.char="#")

# Subsetting for (MuTect's default) HQ filters
somatic_SNVs <- subset(somatic_SNVs, 
                       alt_allele_seen=="True" & 
                         alt_allele_in_normal=="PASS" & 
                         clustered_events=="PASS" &
                         germline_risk == "PASS" & 
                         homologous_mapping_event=="PASS" &
                         multi_event_alt_allele_in_normal=="PASS" &
                         panel_of_normals=="PASS" & 
                         str_contraction=="PASS" &
                         t_lod_fstar=="PASS" &
                         triallelic_site=="PASS" &
                         short_tandem_repeat_membership=="False")

# Reject variants with "Unknown" Hugo Symbols
somatic_SNVs <- subset(somatic_SNVs, Hugo_Symbol != "Unknown")
# Keep only coding mutations
somatic_SNVs <- somatic_SNVs[!somatic_SNVs$Variant_Classification %in% c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", 
                                                                         "IGR", "Intron", "lincRNA", "RNA"),]

genes_df <- data.frame(Gene = unique(somatic_SNVs$Hugo_Symbol),
                                Value = 1, stringsAsFactors = F)

# Exract coding genes with CNV --------------------------------------------
cnv_df <- read.csv("NOTATES/SCNA/all_genes.csv", stringsAsFactors = F)
cnv_df <- cnv_df[cnv_df$ratio <= 0.5 | cnv_df$ratio >= 1.5,]

if(any(duplicated(cnv_df$Gene)))
{
  dups <- names(table(cnv_df$Gene)[table(cnv_df$Gene) > 1])
  
  for(dup in dups)
    cnv_df$ratio[cnv_df$Gene == dup] <- mean(cnv_df$ratio[cnv_df$Gene == dup])
}
cnv_df <- cnv_df[!duplicated(cnv_df$Gene),]

cnv_df$logR <- log2(cnv_df$ratio)
cnv_df <- cnv_df[,c("Gene", "logR")]

final_df <- merge(genes_df, cnv_df, by = "Gene", all = T)
final_df$Value[is.na(final_df$Value)] <- 0
final_df$logR[is.na(final_df$logR)] <- 0

genes_of_interest <- unique(final_df$Gene)

# Enrichment -----------------------------------
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# created named list, eg:  path:map00010: "Glycolysis / Gluconeogenesis"
pathways_list <- keggList("pathway", "hsa")

# make them into KEGG-style human pathway identifiers
pathway_codes <- sub("path:", "", names(pathways_list))

# subsetting by c(TRUE, FALSE) -- which repeats
# as many times as needed, sorts through some
# unexpected packaging of geneIDs in the GENE element
# of each pw[[n]]
genes_by_pathway <- sapply(pathway_codes, function(pwid){
  pw <- keggGet(pwid)
  pw <- pw[[1]]$GENE[c(F, T)]
  pw <- sub(";.+", "", pw)
  pw <- pw[grepl("^[a-zA-Z0-9_-]*$", pw)] ## removing mistaken lines
  pw
})

all_genes <- keys(org.Hs.eg.db, "SYMBOL")

hyperg <- Category:::.doHyperGInternal
hyperg_test <- function(pw_genes, chosen_genes, all_genes, over=TRUE)
{
  pw_genes_selected <- length(intersect(chosen_genes, pw_genes))
  pw_genes_in_pool <- length(pw_genes)
  tol_genes_n_pool <- length(all_genes)
  non_pw_genes_in_pool <- tol_genes_n_pool - pw_genes_in_pool
  num_selected_genes <- length(chosen_genes)
  hyperg(pw_genes_in_pool, non_pw_genes_in_pool,
         num_selected_genes, pw_genes_selected, over)
}

enrichment_res <- t(sapply(genes_by_pathway, hyperg_test, genes_of_interest, all_genes))
enrichment_res <- as.data.frame(enrichment_res)
enrichment_res$p <- unlist(enrichment_res$p)
enrichment_res$odds <- unlist(enrichment_res$odds)
enrichment_res$expected <- unlist(enrichment_res$expected)
enrichment_res <- enrichment_res[order(enrichment_res$p),]
enrichment_res$adj_p <- p.adjust(enrichment_res$p, method = "bonferroni")

# sum(enrichment_res$adj_p < 0.05)
enrichment_res <- enrichment_res[enrichment_res$adj_p < 0.05,]

enrichment_res$Pathway <- pathways_list[match(paste0("path:",rownames(enrichment_res)), names(pathways_list))]
enrichment_res <- enrichment_res[,c(5,1:4)]
enrichment_res$Pathway <- sub(" - Homo sapiens \\(human\\)", "", enrichment_res$Pathway)
enrichment_res <- enrichment_res[order(enrichment_res$adj_p),]

# Visualization -----------------------------------------------------------
suppressPackageStartupMessages(library(pathview))

rownames(final_df) <- final_df$Gene
final_df <- final_df[,-1]
colnames(final_df) <- c("Mutation", "CNV")

dir.create("pathways")
setwd("pathways/")

### Pathways in Cancer
pathview(gene.data = final_df, gene.idtype = "SYMBOL",
         pathway.id = "hsa05200", species = "hsa", out.suffix = "CANCER", 
         keys.align = "y", kegg.native = T, key.pos = "topright", same.layer = F, silent = T)

### Enriched Pathways 
enrichment_res$Somatic_Mutation <- ""
enrichment_res$SCNA_down <- ""
enrichment_res$SCNA_up <- ""

cnv_up <- cnv_df$Gene[cnv_df$logR > 0 ]
cnv_down <- cnv_df$Gene[cnv_df$logR < 0 ]

for(i in 1:nrow(enrichment_res))
{
  tmp <- pathview(gene.data = final_df, gene.idtype = "SYMBOL",
                  pathway.id = rownames(enrichment_res)[i], species = "hsa", out.suffix = enrichment_res$Pathway[i], 
                  keys.align = "y", kegg.native = T, key.pos = "topright", same.layer = F)
  
  tmp <- tmp$plot.data.gene$all.mapped
  tmp <- tmp[tmp!=""]
  tmp <- unlist(strsplit(tmp, split = ","))
  tmp <- unique(tmp)
  tmp <- select(org.Hs.eg.db, tmp, "SYMBOL", "ENTREZID")[,2]
  
  enrichment_res$Somatic_Mutation[i] <- paste(tmp[tmp %in% genes_df$Gene], collapse = ", ")
  enrichment_res$SCNA_down[i] <- paste(tmp[tmp %in% cnv_down], collapse = ", ")
  enrichment_res$SCNA_up[i] <- paste(tmp[tmp %in% cnv_up], collapse = ", ")
}

write.csv(enrichment_res, "enrichment_results.csv", row.names = T)

setwd("..")