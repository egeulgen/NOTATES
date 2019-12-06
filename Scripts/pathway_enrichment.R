##################################################
## Project: NOTATES
## Script purpose: Script for performing KEGG pathway
## enrichment analysis (ORA) using high-confidence genes
## Date: Nov 11, 2019
## Author: Ege Ulgen
##################################################

# Install package(s) if necessary -----------------------------------------
if(!suppressPackageStartupMessages(require(pathfindR))) {
  if (!suppressPackageStartupMessages(require(pak))) {
    install.packages("pak")
  }
  pak::pkg_install("egeulgen/pathfindR")
  suppressPackageStartupMessages(library(pathfindR))
}

options(stringsAsFactors = FALSE)

# High impact somatic SNV/indels ------------------------------------------
somatic_vars <- read.delim("./Oncotator/annotated.sSNVs.tsv", stringsAsFactors=F, comment.char="#")

# Subsetting for (MuTect's default) HQ filters
somatic_vars <- subset(somatic_vars, 
                       alt_allele_seen=="True")

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
  
HQ_mut_res <- run_pathfindR(HQ_mut_df,
                            plot_enrichment_chart = FALSE,
                            visualize_enriched_terms = FALSE,
                            output_dir = "pathfindR_results/HQ_mutations")
write.csv(HQ_mut_res, "pathfindR_results/HQ_mutations/mut_enrichment_results.csv")

### Plot enrichment chart
if (nrow(HQ_mut_res) != 0) {
  HQ_mut_res$Up_regulated <- as.character(HQ_mut_res$Up_regulated)
  HQ_mut_res$Down_regulated <- as.character(HQ_mut_res$Down_regulated)
  png("pathfindR_results/HQ_mutations/enrichment_chart.png", width = 500, height = 400)
  
  if (nrow(HQ_mut_res) > 1) {
    clustered <- cluster_enriched_terms(HQ_mut_res, 
                                        plot_clusters_graph = FALSE, 
                                        plot_dend = FALSE)
    clustered$Up_regulated <- as.character(clustered$Up_regulated)
    clustered$Down_regulated <- as.character(clustered$Down_regulated)
    g <- enrichment_chart(clustered, plot_by_cluster = TRUE)
  } else {
    g <- enrichment_chart(HQ_mut_res, plot_by_cluster = TRUE)
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
                         PVAL = 0.05)

HQ_SCNA_res <- run_pathfindR(HQ_SCNA_df,
                            plot_enrichment_chart = FALSE,
                            visualize_enriched_terms = FALSE,
                            output_dir = "pathfindR_results/HQ_SCNA")
write.csv(HQ_SCNA_res, "pathfindR_results/HQ_SCNA/SCNA_enrichment_results.csv")

### Plot enrichment chart\
if (nrow(HQ_SCNA_res) != 0) {
  HQ_SCNA_res$Up_regulated <- as.character(HQ_SCNA_res$Up_regulated)
  HQ_SCNA_res$Down_regulated <- as.character(HQ_SCNA_res$Down_regulated)
  png("./pathfindR_results/HQ_SCNA/enrichment_chart.png", width = 500, height = 700)
  if (nrow(HQ_SCNA_res) > 1) {
    clustered <- cluster_enriched_terms(HQ_SCNA_res, 
                                        plot_clusters_graph = FALSE, 
                                        plot_dend = FALSE)
    g <- enrichment_chart(clustered, plot_by_cluster = TRUE)
  } else {
    g <- enrichment_chart(HQ_SCNA_res, plot_by_cluster = TRUE)
  }
  print(g)
  dev.off()
}

