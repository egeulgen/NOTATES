############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                    Germline Analysis                     #
#                   Ege Ulgen, Dec 2017                    #
############################################################

# Locate directory for data sources Set working dir and create output dir ----------------------

# dir for data sources (same as the directory where script is located)
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script_dir <- dirname(script.name)
  
# set working directory
setwd("./Germline")
dir.create("output")

# 0. Seperate Report for Common variants with low penetrance --------------
# load annotated SNV file
germline <- read.delim("../Oncotator/annotated.germline_SNVs.tsv", comment.char = "#", stringsAsFactors = F)

# discard if alt. allele not seen
germline <- subset(germline, alt_allele_seen == "True")

# Subsetting for Germline HQ filters
germline <- subset(germline, HQ_SNP_filter=="PASS" & HQ_InDel_filter=="PASS" & LowQual=="PASS")

# load gCCV table
common_vars <- read.csv(paste0(script_dir, "/common_variants_list.csv"), stringsAsFactors = F)

# Create df for variants that predispose to glioma in the sample
idx <- which(germline$id %in% common_vars$rs_ID)
common_variant_df <- germline[idx, ]
# make sure alt. allele was seen
common_variant_df <- common_variant_df[common_variant_df$alt_allele_seen == "True", ]

cols <- c("id", "Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification",
          "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "read_depth", "allelic_depth", 
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
rm(list=setdiff(ls(), c("germline","script_dir")))

# Preprocess --------------------------------------------------------------
# Further subsetting 
germline <- subset(germline, read_depth>20) #### any other?????

# # Filter out non-coding variants
# rejects <- c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", 
#              "Intron","lincRNA", "RNA", "De_novo_Start_InFrame")
# 
# germline <- subset(germline, !(Variant_Classification %in% rejects))

# I. Extract relevant genes with preset filter groups ---------------------
filter_df <- read.csv(paste0(script_dir, "/Germline_filter_list.csv"), stringsAsFactors = F)

cgc_df <- read.csv(paste0(script_dir,"/../CGC_dec15_17.csv"), stringsAsFactors = F)
cgc_df <- data.frame(Gene.Name = cgc_df$Gene.Symbol,
                     Comment = "",
                     Group = "Cancer Gene Census",
                     Source = "COSMIC",
                     Rank = 3,
                     short_name = "CGC", stringsAsFactors = F)

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

write.csv(germline_final,"./output/germline_variant_NO_RS_FILTER.csv", row.names = F)

# II. Evaluate Variants ---------------------------------------------------
pathogenic_SNPs <- read.table(paste0(script_dir, "/dbSNP_table.txt"), stringsAsFactors = F, header = T)

keep <- germline_final$id %in% pathogenic_SNPs$rs_id

germline_final <- germline_final[keep, ]

# III. Report Relevant Variants -------------------------------------------
cols_to_keep <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", 
                  "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "id", 
                  "Filter_Group", "Filter_Comment",
                  "allelic_depth", "allele_frequency")

germline_final <- germline_final[, cols_to_keep]
colnames(germline_final) <- gsub("Tumor", "Germline", colnames(germline_final))
germline_final$minor_allele <- pathogenic_SNPs$Minor_allele[match(germline_final$id,pathogenic_SNPs$rs_id)]

write.csv(germline_final,"./output/germline_variant_report.csv", row.names = F)
