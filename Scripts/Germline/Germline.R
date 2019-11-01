##################################################
## Project: NOTATES
## Script purpose: Script for processing germline
## SNP/indels for categorical reporting of variants
## that are reported as "pathogenic"/"likely-pathogenic"
## in ClinVar
## Date: Oct 23, 2019
## Author: Ege Ulgen
##################################################

# Locate directory for data sources Set working dir and create output dir ----------------------

# dir for data sources (same as the directory where script is located)
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script_dir <- dirname(script.name)

options(stringsAsFactors = FALSE)
  
# set working directory
setwd("./Germline")
dir.create("output")

# Load oncotator-annotated germline SNPs ----------------------------------
# load annotated SNV file
germline <- read.delim("../Oncotator/annotated.germline_SNVs.tsv", comment.char = "#")

# discard if alt. allele not seen
germline <- subset(germline, alt_allele_seen == "True")

# Subsetting for Germline HQ filters
germline <- subset(germline, HQ_SNP_filter=="PASS" & HQ_InDel_filter=="PASS" & LowQual=="PASS")

# Add ref/alt depths
germline$Ref_depth <- vapply(germline$allelic_depth, function(x) 
  as.integer(unlist(strsplit(x, split =","))[1]), 1L)
germline$Alt_depth <- vapply(germline$allelic_depth, function(x) 
  as.integer(unlist(strsplit(x, split =","))[2]), 1L)  
  
# 0. Seperate Report for Common variants with low penetrance --------------
# load gCCV table
common_vars <- read.csv(paste0(script_dir, "/common_variants_list.csv"))

# Create df for variants that predispose to glioma in the sample
idx <- which(germline$id %in% common_vars$rs_ID)
common_variant_df <- germline[idx, ]

cols <- c("id", "Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification",
          "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Ref_depth", "Alt_depth", 
          "allele_frequency")

common_variant_df <- common_variant_df[, cols]
colnames(common_variant_df) <- gsub("Tumor", "Normal", colnames(common_variant_df))

# add further information
idx <- match(common_variant_df$id, common_vars$rs_ID)
common_variant_df$gCCV_functional_consequence <- common_vars$Functional_consequence[idx]
common_variant_df$gCCV_Risk_allele <- common_vars$Risk_allele[idx]
common_variant_df$gCCV_Ancestral_allele <- common_vars$Ancestral_allele[idx]
common_variant_df$gCCV_MAF_1000G <- common_vars$MAF_1kG[idx]

write.csv(common_variant_df, "./output/common_variant_report.csv", row.names = F)

# Preprocess --------------------------------------------------------------
# I. Extract relevant genes with preset filter groups ---------------------
filter_df <- read.csv(paste0(script_dir, "/Germline_filter_list.csv"))

cgc_df <- read.csv(paste0(script_dir,"/../CGC_latest.csv"))
cgc_df <- data.frame(Gene.Name = cgc_df$Gene.Symbol,
                     Comment = "",
                     Group = "Cancer Gene Census",
                     Source = "COSMIC",
                     Rank = 3,
                     short_name = "CGC")

filter_df <- rbind(filter_df, cgc_df)

# lookup indices of genes/variants to be filtered
idx <- which(germline$Hugo_Symbol %in% filter_df$Gene.Name)
germline_final <- germline[idx, ]

# annotate filter groups
germline_final$Rank <- 0
germline_final$Filter_Group <- NA
germline_final$Filter_Comment <- NA

for (i in 1:nrow(germline_final)) {
  tmp_gname <- germline_final$Hugo_Symbol[i]
  tmp <- filter_df[filter_df$Gene.Name == tmp_gname, ]
  # rank - highest
  germline_final$Rank[i] <- max(tmp$Rank)
  # filter group(s)
  germline_final$Filter_Group[i] <- paste(tmp$short_name, collapse = ", ")
  # any comment(s)
  tmp_comment <- tmp$Comment
  tmp_comment <- tmp_comment[tmp_comment != ""]
  germline_final$Filter_Comment[i] <- paste(tmp_comment, collapse = ", ")
}
germline_final <- germline_final[order(germline_final$Rank), ]
rm(list = setdiff(ls(), c("script_dir", "germline_final")))

write.csv(germline_final,"./output/germline_variants_NO_FILTER.csv", row.names = F)

# II. Evaluate Variants ---------------------------------------------------
clinvar_pathogenic <- read.table(file.path(script_dir, 
                                           "ClinVar_pathogenic_positions.txt"), header = FALSE)
# both 1-based
germline_final$lookup <- paste(germline_final$Chromosome, 
                               germline_final$Start_position, 
                               germline_final$Reference_Allele, 
                               germline_final$Tumor_Seq_Allele2, 
                               sep = "_")
clinvar_pathogenic$lookup <- paste(clinvar_pathogenic$V1,
                                   clinvar_pathogenic$V2,
                                   clinvar_pathogenic$V4,
                                   clinvar_pathogenic$V5,
                                   sep = "_")

germline_final <- germline_final[germline_final$lookup %in% clinvar_pathogenic$lookup, ]

# III. Report Relevant Variants -------------------------------------------
cols_to_keep <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", 
                  "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "id", 
                  "Filter_Group", "Filter_Comment",
                  "Ref_depth", "Alt_depth", "allele_frequency")

germline_final <- germline_final[, cols_to_keep]
colnames(germline_final) <- gsub("Tumor", "Germline", colnames(germline_final))

write.csv(germline_final,"./output/germline_variant_report.csv", row.names = F)
