##################################################
## Project: NOTATES
## Script purpose: Script for performing KEGG pathway
## enrichment analysis with pathfindR using high-confidence genes
## Date: Sep 23, 2020
## Author: Ege Ulgen
##################################################

# Install package(s) if necessary -----------------------------------------
if(!suppressPackageStartupMessages(require(pathfindR))) {
  devtools::install_github("egeulgen/pathfindR")
  suppressPackageStartupMessages(library(pathfindR))
}

# High impact somatic SNV/indels ------------------------------------------
somatic_vars <- read.delim("./Funcotator/annotated_somatic.maf", stringsAsFactors=F, comment.char="#")

# filter for VAF > 0.05
somatic_vars <- somatic_vars[somatic_vars$tumor_f > 0.05, ]

# Exclude variants in FLAGs
flags <- c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", 
           "MUC5B", "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", 
           "LAMA5", "AHNAK", "HMCN1", "USH2A", "DNAH11", "MACF1", 
           "MUC17")
somatic_vars <- somatic_vars[!somatic_vars$Hugo_Symbol %in% flags, ]

# Only include somatic snv/indels with Variant Classifications of High/Moderate variant consequences
high_conseq <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                 "Translation_Start_Site","Nonsense_Mutation", 
                 "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", 
                 "Missense_Mutation")
somatic_vars <- somatic_vars[somatic_vars$Variant_Classification %in% high_conseq, ]

HQ_mut_df <- data.frame(GENE = somatic_vars$Hugo_Symbol,
                        PVAL = 0.05)

HQ_mut_res <- tryCatch({
  res <- run_pathfindR(HQ_mut_df,
                       plot_enrichment_chart = FALSE,
                       visualize_enriched_terms = FALSE,
                       output_dir = "pathfindR_results/HQ_mutations")
  res
}, error = function(e) {
  message("Cannot perform SNV/indel enrichment analysis!")
  message("Here's the original error message:")
  message(e)
  return(data.frame())
})
write.csv(HQ_mut_res, "pathfindR_results/HQ_mutations/mut_enrichment_results.csv")

### Plot enrichment chart
if (nrow(HQ_mut_res) != 0) {
  HQ_mut_res$Up_regulated <- as.character(HQ_mut_res$Up_regulated)
  HQ_mut_res$Down_regulated <- as.character(HQ_mut_res$Down_regulated)
  png("pathfindR_results/HQ_mutations/enrichment_chart.png", width = 500, height = 400)
  
  if (nrow(HQ_mut_res) > 2) {
    clustered <- cluster_enriched_terms(HQ_mut_res, 
                                        plot_clusters_graph = FALSE, 
                                        plot_dend = FALSE)
    clustered$Up_regulated <- as.character(clustered$Up_regulated)
    clustered$Down_regulated <- as.character(clustered$Down_regulated)
    g <- enrichment_chart(clustered, plot_by_cluster = TRUE)
  } else {
    g <- enrichment_chart(HQ_mut_res)
  }
  print(g)
  dev.off()
}


# High impact SCNAs -------------------------------------------------------
cnv_df <- read.csv("NOTATES/SCNA/all_genes.csv", stringsAsFactors = F)
cnv_df <- cnv_df[cnv_df$ratio <= 0.5 | cnv_df$ratio >= 1.5, ]

if (any(duplicated(cnv_df$Gene))) {
  tmp_df <- c()
  for (i in seq_len(nrow(cnv_df))) {
    if (sum(cnv_df$Gene == cnv_df$Gene[i]) == 1) {
      tmp_df <- rbind(tmp_df, cnv_df[i, ])
    } else {
      tmp <- cnv_df[cnv_df$Gene == cnv_df$Gene[i], ]
      tmp$ratio <- mean(tmp$ratio)
      tmp_df <- rbind(tmp_df, tmp[1, ])
    }
  }
  cnv_df <- tmp_df
}

HQ_SCNA_df <- data.frame(Gene = cnv_df$Gene,
                         Change = ifelse(cnv_df$ratio > 1, 1, -1),
                         PVAL = rep(0.05, length(cnv_df$Gene)))

HQ_SCNA_res <- tryCatch({
  res <- run_pathfindR(HQ_SCNA_df,
                       plot_enrichment_chart = FALSE,
                       visualize_enriched_terms = FALSE,
                       output_dir = "pathfindR_results/HQ_SCNA")
  res
}, error = function(e) {
  message("Cannot perform SCNA enrichment analysis!")
  message("Here's the original error message:")
  message(e)
  return(data.frame())
})
write.csv(HQ_SCNA_res, "pathfindR_results/HQ_SCNA/SCNA_enrichment_results.csv")

### Plot enrichment chart\
if (nrow(HQ_SCNA_res) != 0) {
  HQ_SCNA_res$Up_regulated <- as.character(HQ_SCNA_res$Up_regulated)
  HQ_SCNA_res$Down_regulated <- as.character(HQ_SCNA_res$Down_regulated)
  png("./pathfindR_results/HQ_SCNA/enrichment_chart.png", width = 500, height = 700)
  if (nrow(HQ_SCNA_res) > 2) {
    clustered <- cluster_enriched_terms(HQ_SCNA_res, 
                                        plot_clusters_graph = FALSE, 
                                        plot_dend = FALSE)
    g <- enrichment_chart(clustered, plot_by_cluster = TRUE)
  } else {
    g <- enrichment_chart(HQ_SCNA_res)
  }
  print(g)
  dev.off()
}
