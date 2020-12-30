##################################################
## Project: NOTATES
## Script purpose: Script for processing germline
## SNP/indels for categorical reporting of variants
## that:
## 1. are NOT reported as "benign"/"likely-benign"
## in ClinVar
## 2. have MAF < 1%
## 3. are non-synonymous
## 4. are not in FLAGs
## Date: Dec 7, 2020
## Author: Ege Ulgen
##################################################

# Locate directory for data sources Set working dir and create output dir ----------------------

# dir for data sources (same as the directory where script is located)
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script_dir <- dirname(script.name)

# set working directory
setwd("./Germline")
dir.create("output")

# Load annotated germline variants ----------------------------------------
# load annotated SNV file
germline <- read.delim(file.path(dirname(getwd()), "Funcotator", 
                                 "annotated_germline.maf"), 
                       comment.char = "#")

# Add ref/alt depths
germline$Ref_depth <- germline$t_ref_count
germline$Alt_depth <- germline$t_alt_count  

# add id column
germline$id <- ifelse(germline$gnomAD_exome_ID != "", germline$gnomAD_exome_ID, germline$gnomAD_genome_ID)

germline$Clin_sig <- germline$ClinVar_VCF_CLNSIG
germline$Clin_sig[is.na(germline$Clin_sig)] <- ""

# fix long clin. sig.s
sig_desc <- "Conflicting_interpretations_of_pathogenicity"
germline$Clin_sig[grepl(sig_desc, germline$Clin_sig)] <- gsub(sig_desc, "Conflicting", germline$Clin_sig[grepl(sig_desc, germline$Clin_sig)] )

sig_desc <- "Uncertain_significance"
germline$Clin_sig[grepl(sig_desc, germline$Clin_sig)] <- gsub(sig_desc, "VUS", germline$Clin_sig[grepl(sig_desc, germline$Clin_sig)] )

sig_desc <- "not_provided"
germline$Clin_sig[grepl(sig_desc, germline$Clin_sig)] <- ""

germline$Clin_sig[germline$Clin_sig == ""] <- "not reported" 

# 0. Seperate Report for Common variants with low penetrance --------------
# load gCCV table
common_vars <- read.csv(paste0(script_dir, "/common_variants_list.csv"))

# Create df for variants that predispose to glioma in the sample
common_variant_df <- germline[germline$id %in% common_vars$rs_ID, ]

cols <- c("id", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification",
          "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Ref_depth", "Alt_depth", 
          "AF")

common_variant_df <- common_variant_df[, cols]
colnames(common_variant_df) <- gsub("Tumor", "Normal", colnames(common_variant_df))

# add further information
idx <- match(common_variant_df$id, common_vars$rs_ID)
common_variant_df$gCCV_functional_consequence <- common_vars$Functional_consequence[idx]
common_variant_df$gCCV_Risk_allele <- common_vars$Risk_allele[idx]
common_variant_df$gCCV_Ancestral_allele <- common_vars$Ancestral_allele[idx]
common_variant_df$gCCV_MAF_1000G <- common_vars$MAF_1kG[idx]

write.csv(common_variant_df, "output/common_variant_report.csv", row.names = F)

# I. Extract relevant genes with preset filter groups ---------------------
filter_df <- read.csv(file.path(script_dir, "Germline_filter_list.csv"))

main_scripts_dir <- dirname(script_dir)

cgc_df <- read.csv(file.path(main_scripts_dir, "CGC_latest.csv"))
cgc_df <- data.frame(Gene.Name = cgc_df$Gene.Symbol,
                     Comment = "",
                     Group = "Cancer Gene Census",
                     Source = "COSMIC",
                     Rank = 3,
                     short_name = "CGC")
filter_df <- rbind(filter_df, cgc_df)

ddr_df <- read.csv(file.path(main_scripts_dir, "NOTATES", "curated_dbs", "DNA_damage_repair_1Jul20.csv"))
ddr_df <- data.frame(Gene.Name = ddr_df$Gene.Name,
                     Comment = ddr_df$FUNCTION,
                     Group = "DNA Damage Repair Gene",
                     Source = "https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html",
                     Rank = 2,
                     short_name = "DDR")
filter_df <- rbind(filter_df, ddr_df)

# lookup indices of genes/variants to be filtered
germline_final <- germline[germline$Hugo_Symbol %in% filter_df$Gene.Name, ]

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

write.csv(germline_final,"./output/germline_variants_NO_FILTER.csv", row.names = F)

# II. Evaluate Variants ---------------------------------------------------
# only include pathogenic variants
germline_final <- germline_final[grepl("pathogenic", germline_final$Clin_sig, ignore.case = TRUE), ]

# Exclude variants in FLAGs
flags <- c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG",
           "MUC5B", "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2",
           "LAMA5", "AHNAK", "HMCN1", "USH2A", "DNAH11", "MACF1",
           "MUC17")
germline_final <- germline_final[!germline_final$Hugo_Symbol %in% flags, ]

# III. Report Relevant Variants -------------------------------------------
cols_to_keep <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", 
                  "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "id", 
                  "Filter_Group", "Filter_Comment",
                  "Ref_depth", "Alt_depth", "AF", "Clin_sig")

germline_final <- germline_final[, cols_to_keep]
colnames(germline_final) <- gsub("Tumor", "Germline", colnames(germline_final))

write.csv(germline_final,"./output/germline_variant_report.csv", row.names = FALSE)
