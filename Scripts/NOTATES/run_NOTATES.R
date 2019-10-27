##################################################
## Project: NOTATES
## Script purpose: Script for sequentially filtering
## germline and somatic alterations for reporting of
## clinically-relevant findings
## Date: Oct 23, 2019
## Author: Ege Ulgen
##################################################

dir.create("./NOTATES/")
setwd("./NOTATES/")

options(stringsAsFactors = FALSE)

# dir for data sources (same as the directory where script is located)
initial_options <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", initial_options[grep("--file=", initial_options)])
notates_dir <- dirname(script_name)

## Argument for MB vs Glioma
mb_arg <- commandArgs(trailingOnly=TRUE)
if (length(mb_arg) == 0)
  mb_arg[1] = 'glioma'

# Necessary resources -----------------------------------------------------
# CGC
CGC_df <- read.csv(paste0(notates_dir, "/../CGC_latest.csv"))

# Curated Databases
dna_repair_df <- read.csv(paste0(notates_dir, "/curated_dbs/DNA_damage_repair_22feb16.csv"))
irinotecan_df <- read.csv(paste0(notates_dir, "/curated_dbs/irinotecan_response_genes_22feb16.csv"))
tmz_df <- read.csv(paste0(notates_dir, "/curated_dbs/temozolamide_resistance_genes_22feb16.csv"))
if (mb_arg == 'glioma') {
  KEGG_df <- read.csv(paste0(notates_dir, "/curated_dbs/important_KEGG_pws_19mar19.csv"))
} else {
  KEGG_df <- read.csv(paste0(notates_dir, "/curated_dbs/important_KEGG_pws_MB_19mar19.csv"))
}

# Curated Alterations
if (mb_arg == 'glioma') {
  curated_SNV <- read.csv(paste0(notates_dir, "/curated_alterations/glioma_important_SNV_May25_17.csv"))
  curated_CNA <- read.csv(paste0(notates_dir, "/curated_alterations/glioma_important_CNA_May25_17.csv"))
} else {
  curated_SNV <- read.csv(paste0(notates_dir, "/curated_alterations/MB_important_SNV_Mar19_19.csv"))
  curated_SNV$Genomic_Alt[is.na(curated_SNV$Genomic_Alt)] <- ""
  curated_SNV$Protein_Alt[is.na(curated_SNV$Protein_Alt)] <- ""
  # keep only >= 5 times affected
  curated_SNV <- curated_SNV[curated_SNV$Note >= 5, ]
  curated_CNA <- read.csv(paste0(notates_dir, "/curated_alterations/MB_important_CNA_Mar19_19.csv"))
}

# Germline Mutations ------------------------------------------------------
dir.create("./Germline")
germline_mutations <- read.csv("../Germline/output/germline_variant_report.csv")

germline_cols <- c("Hugo_Symbol", "id", "Chromosome", "Start_position",
                   "Filter_Comment", "Variant_Classification",
                   "Reference_Allele", "Germline_Seq_Allele2", "allelic_depth")
renamed_g_cols <- c("Gene", "rs_id", "Chr", "Pos", 
                    "Disease(s)", "Effect", 
                    "Ref", "Alt", "AD")

### ACMG incidental
if(any(grepl("ACMG",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("ACMG", germline_mutations$Filter_Group),]
  tmp <- tmp[, germline_cols]
  colnames(tmp) <- renamed_g_cols
  germline_mutations <- germline_mutations[!grepl("ACMG", germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/ACMG.csv", row.names = F)
}
### Cancer Gene Census
if(any(grepl("CGC",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("CGC",germline_mutations$Filter_Group),]
  tmp <- tmp[, germline_cols]
  colnames(tmp) <- renamed_g_cols
  tmp <- tmp[, colnames(tmp) != "Disease(s)"]
  
  germline_mutations <- germline_mutations[!grepl("CGC",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/CGC.csv", row.names = F)
}
### Cancer Predisposition Gene
if(any(grepl("CPG",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("CPG",germline_mutations$Filter_Group),]
  tmp <- tmp[, germline_cols]
  colnames(tmp) <- renamed_g_cols
  tmp <- tmp[, colnames(tmp) != "Disease(s)"]
  germline_mutations <- germline_mutations[!grepl("CPG",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/CPG.csv", row.names = F)
}
### Fanconi Anemia Pathway
if(any(grepl("FAP",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("FAP",germline_mutations$Filter_Group),]
  tmp <- tmp[, germline_cols]
  colnames(tmp) <- renamed_g_cols
  tmp <- tmp[, colnames(tmp) != "Disease(s)"]
  germline_mutations <- germline_mutations[!grepl("FAP",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/FAP.csv", row.names = F)
}
### Other
if(nrow(germline_mutations) != 0)
{
  tmp <- germline_mutations
  tmp <- tmp[, germline_cols]
  colnames(tmp) <- renamed_g_cols
  tmp <- tmp[, colnames(tmp) != "Disease(s)"]
  
  write.csv(tmp, "./Germline/OTHER.csv", row.names = FALSE)
}

### Common Variants
common_var <- read.csv("../Germline/output/common_variant_report.csv")
common_var <- common_var[,c("id","Hugo_Symbol","Variant_Classification",
                            "gCCV_Risk_allele","Reference_Allele","Normal_Seq_Allele2","allelic_depth")]
colnames(common_var) <-  c("rs_id","Gene","Effect","Risk Allele", "Ref", "Alt", "AD")
write.csv(common_var, "./Germline/common_var.csv", row.names = FALSE)

# Somatic Mutations -------------------------------------------------------
dir.create("./Somatic_SNV")
somatic_SNVs <- read.delim("../Oncotator/annotated.sSNVs.tsv", comment.char="#")
somatic_SNVs <- somatic_SNVs[order(somatic_SNVs$tumor_f, decreasing = TRUE), ]

# Subsetting for (MuTect's default) HQ filters
somatic_SNVs <- subset(somatic_SNVs, 
                       alt_allele_seen=="True" & 
                         short_tandem_repeat_membership == "False")

# # Reject variants with "Unknown" Hugo Symbols
# somatic_SNVs <- subset(somatic_SNVs, Hugo_Symbol != "Unknown")

# Change "" protein changes to NA
somatic_SNVs$Protein_Change[somatic_SNVs$Protein_Change == ""] <- NA

# keep only needed annotations
keep <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", "Variant_Type", 
          "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "tumor_f", "genotype", "dbSNP_RS", "Genome_Change", "Codon_Change",
          "Protein_Change", "UniProt_AApos", "DNARepairGenes_Role", "FamilialCancerDatabase_Syndromes",
          "COSMIC_n_overlapping_mutations", "COSMIC_total_alterations_in_gene", "dbNSFP_SIFT_pred", 
          colnames(somatic_SNVs)[grepl("^GO_", colnames(somatic_SNVs))], 
          colnames(somatic_SNVs)[grepl("^CGC_", colnames(somatic_SNVs))])

write.csv(somatic_SNVs[, keep], "./Somatic_SNV/detailed_sSNV.csv", row.names = F, quote = F)

somatic_SNVs <- somatic_SNVs[, c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "Genome_Change",
                                 "tumor_f", "genotype", "COSMIC_n_overlapping_mutations",
                                 "COSMIC_total_alterations_in_gene","UniProt_Region","dbNSFP_SIFT_pred")]

noncoding_SNVs <- somatic_SNVs[somatic_SNVs$Variant_Classification %in% c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", 
                                                                          "IGR", "Intron", "lincRNA", "RNA"),]
noncoding_SNVs <- noncoding_SNVs[,setdiff(colnames(noncoding_SNVs),c("Protein_Change","UniProt_Region"))]
write.csv(noncoding_SNVs, "Somatic_SNV/noncoding_SNVs.csv", row.names = F)

somatic_SNVs <- somatic_SNVs[!somatic_SNVs$Variant_Classification %in% c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", 
                                                                         "IGR", "Intron", "lincRNA", "RNA"),]
#### Other annotations 
somatic_SNVs$DNA_repair <- ifelse(somatic_SNVs$Hugo_Symbol %in% dna_repair_df$Gene.Name, "yes", "no")
somatic_SNVs$DNA_repair[somatic_SNVs$DNA_repair == "yes"] <- dna_repair_df$FUNCTION[match(somatic_SNVs$Hugo_Symbol[somatic_SNVs$DNA_repair == "yes"],dna_repair_df$Gene.Name)]

somatic_SNVs$Irinotecan_resp <- ifelse(somatic_SNVs$Hugo_Symbol %in% irinotecan_df$Genes, "yes", "no")
somatic_SNVs$TMZ_resistance <- ifelse(somatic_SNVs$Hugo_Symbol %in% tmz_df$gene, "yes", "no")

somatic_SNVs$selected_KEGG <- ifelse(somatic_SNVs$Hugo_Symbol %in% KEGG_df$Gene, "yes", "no")
somatic_SNVs$selected_KEGG[somatic_SNVs$selected_KEGG == "yes"] <- KEGG_df$pathway[match(somatic_SNVs$Hugo_Symbol[somatic_SNVs$selected_KEGG == "yes"],KEGG_df$Gene)]


### Important Glioma SNVs
if(any(somatic_SNVs$Hugo_Symbol %in% curated_SNV$Gene))
{
  idx <- which(somatic_SNVs$Hugo_Symbol %in% curated_SNV$Gene)
  keep <- c()
  for(i in idx)
  {
    tmp <- somatic_SNVs[i, ]
    idx2 <- which(curated_SNV$Gene == tmp$Hugo_Symbol)
    
    if(any(curated_SNV$Genomic_Alt[idx2] != ""))
    {
      tmp2 <- curated_SNV$Genomic_Alt[idx2]
      if(unlist(strsplit(tmp$Genome_Change,":"))[2] %in% tmp2)
        keep <- c(keep, i)
    }
    else if(any(curated_SNV$Protein_Alt[idx2] != ""))
    {
      tmp2 <- curated_SNV$Protein_Alt[idx2]
      if(unlist(strsplit(tmp$Protein_Change,"\\."))[2] %in% tmp2)
        keep <- c(keep, i)
      else if(any(grepl("\\*",tmp2)))
      {
        tmp2 <- tmp2[grepl("\\*",tmp2)]
        tmp2 <- sapply(tmp2, function(x) substr(x, 1, nchar(x)-1))
        
        tmp3 <- unlist(strsplit(tmp$Protein_Change,"\\."))[2]
        if(substr(tmp3,1,nchar(tmp3)-1) %in% tmp2)
          keep <- c(keep, i)
      }
    }
    else
      keep <- c(keep, i)
  }
  tmp <- somatic_SNVs[keep,]
  somatic_SNVs <- somatic_SNVs[-keep,]
  write.csv(tmp, "Somatic_SNV/important_SNVs.csv",row.names = F)
}

### in CGC and COSMIC hotspot
if(any(somatic_SNVs$Hugo_Symbol %in% CGC_df$Gene.Symbol))
{
  cgc_id <- which(somatic_SNVs$Hugo_Symbol %in% CGC_df$Gene.Symbol)
  hotspot_id <- which(somatic_SNVs$Hugo_Symbol %in% CGC_df$Gene.Symbol & somatic_SNVs$COSMIC_n_overlapping_mutations>0) # should this be >10 or more?
  if(length(hotspot_id)>0)
  {
    cgc_id <- setdiff(cgc_id,hotspot_id)
    tmp <- somatic_SNVs[hotspot_id,]
    write.csv(tmp, "Somatic_SNV/COSMIC_hotspot.csv", row.names = F)
  }
  
  if(length(cgc_id) >0)
  {
    tmp <- somatic_SNVs[cgc_id,]
    write.csv(tmp, "Somatic_SNV/CGC_genes.csv", row.names = F)
  }
  somatic_SNVs <- somatic_SNVs[-c(cgc_id,hotspot_id),]
}

### DNA damage repair
if(any(somatic_SNVs$DNA_repair != "no"))
{
  tmp <- somatic_SNVs[somatic_SNVs$DNA_repair != "no",]
  somatic_SNVs <- somatic_SNVs[somatic_SNVs$DNA_repair == "no",]
  
  write.csv(tmp, "Somatic_SNV/DDR_related.csv",row.names = F)
}

### Important KEGG Pathways
if(any(somatic_SNVs$selected_KEGG != "no"))
{
  tmp <- somatic_SNVs[somatic_SNVs$selected_KEGG != "no",]
  somatic_SNVs <- somatic_SNVs[somatic_SNVs$selected_KEGG == "no",]
  
  write.csv(tmp, "Somatic_SNV/KEGG_selected.csv",row.names = F)
}

### Irinotecan Response
if(any(somatic_SNVs$Irinotecan_resp != "no"))
{
  tmp <- somatic_SNVs[somatic_SNVs$Irinotecan_resp != "no",]
  somatic_SNVs <- somatic_SNVs[somatic_SNVs$Irinotecan_resp == "no",]
  
  write.csv(tmp, "Somatic_SNV/irinotecan.csv",row.names = F)
}

### TMZ Resistance
if(any(somatic_SNVs$TMZ_resistance != "no"))
{
  tmp <- somatic_SNVs[somatic_SNVs$TMZ_resistance != "no",]
  somatic_SNVs <- somatic_SNVs[somatic_SNVs$TMZ_resistance == "no",]
  
  write.csv(tmp, "Somatic_SNV/tmz_res.csv",row.names = F)
}
### Other coding
write.csv(somatic_SNVs, "Somatic_SNV/Other_coding.csv", row.names = F)


# Helper functions for Cytoband & Gene annotations ------------------------
source(file.path(notates_dir, "segment_annotation_hg19.R"))

# read in cytobands file
cytobands_df <- read.delim(file.path(notates_dir, "hg19_cytoBand.txt"), 
                           header = FALSE)
colnames(cytobands_df) <- c("Chr", "Start", "End", "Cytb_name", "stain")
cytobands_df$Cytb_name <- paste0(cytobands_df$Chr, cytobands_df$Cytb_name)

# CNV - exomeCNV ---------------------------------------------------------------------
dir.create("./SCNA")
# read in cnv file
cnv <- read.delim("../ExomeCNV/CNV.cnv.txt", header = TRUE)
# discard rows with NA spec
cnv <- subset(cnv, !is.na(spec))
# discard rows with NA sens
cnv <- subset(cnv, !is.na(sens))

# -0.25 and 0.2 for cut-off
#The inner cutoffs of +0.2 and -0.25 are sensitive 
#enough to detect a single-copy gain or loss in a 
#diploid tumor with purity (or subclone cellularity) as low as 30%. 

cnv <- cnv[cnv$logR <= -0.25 | cnv$logR >= 0.2,]
# cnv <- cnv[cnv$ratio <= 0.5 | cnv$ratio >= 1.5,]

copy_num_call <- function(ratio){
  if( ratio < 0.5 )
    return(0)
  else if( ratio < 1 )
    return(1)
  else
    return(round(ratio*2))
}

cnv$copy.number <- vapply(cnv$ratio, copy_num_call, 2)

## annotate genes and cytobands within segments
cnv_genes_overlap <- annotate_genes(cnv)
cnv_cytb_overlap <- annotate_cytb(cnv, cytobands_df)

cnv$length <- cnv$probe_end - cnv$probe_start + 1
cnv$genes <- sapply(cnv_genes_overlap$Segment_genes, function(x) ifelse(length(x) > 100, ">100", paste(x, collapse = ", ")))
cnv$cytoband <- sapply(cnv_cytb_overlap$Segment_cytbs, function(cyt) shorten_cyt_anno(cyt, cytobands_df))

cnv <- cnv[, c("genes", "length", "cytoband", setdiff(colnames(cnv), c("genes", "cytoband", "length")))]
write.csv(cnv, "SCNA/exomeCNV.csv", row.names = FALSE)

cnv_by_gene <- data.frame(Gene=names(cnv_genes_overlap$Gene_segments),
                          Segment = NA, ratio = NA, CN = NA, av_cov = NA)
cnv$id <- paste0(cnv$chr,":", cnv$probe_start, "-", cnv$probe_end)
gene_segs <- cnv_genes_overlap$Gene_segments

original_len <- length(gene_segs)
for(i in seq_len(original_len)) {
  segments_vec <- gene_segs[[i]]
  idx <- match(segments_vec, cnv$id)
  
  if(length(segments_vec) == 1) {
    cnv_by_gene$ratio[i] <- cnv$ratio[idx]
    cnv_by_gene$av_cov[i] <- cnv$average.coverage[idx]
    cnv_by_gene$Segment[i] <- segments_vec
  } else {
    cnv_by_gene$ratio[i] <- cnv$ratio[idx[1]]
    cnv_by_gene$av_cov[i] <- cnv$average.coverage[idx[1]]
    cnv_by_gene$Segment[i] <- names(gene)[1]
    for(j in 2:length(segments_vec))
    {
      cnv_by_gene <- rbind(cnv_by_gene, 
                           data.frame(Gene = names(gene_segs)[i], 
                                      Segment = segments_vec[j],
                                      ratio = cnv$ratio[idx[j]],
                                      CN = NA,
                                      av_cov = cnv$average.coverage[idx[j]]))
    }
  }
}
cnv_by_gene$av_cov <- as.numeric(cnv_by_gene$av_cov)
cnv_by_gene$ratio <- as.numeric(cnv_by_gene$ratio)
cnv_by_gene$CN <- sapply(cnv_by_gene$ratio, copy_num_call)

write.csv(cnv_by_gene, "SCNA/all_genes.csv", row.names = F)

### Glioma important SCNAs
if(any(cnv_by_gene$Gene %in% curated_CNA$Gene)) {
  tmp <- cnv_by_gene[cnv_by_gene$Gene %in% curated_CNA$Gene,]
  tmp <- tmp[tmp$CN != 2,]
  if (nrow(tmp) != 0) {
    tmp$keep <- T
    for(i in 1:nrow(tmp)) {
      tmp2 <- tmp[i,]
      tmp3 <- curated_CNA$Type[curated_CNA$Gene == tmp2$Gene]
      if( !(tmp3 == "A" & tmp2$ratio >1) & !(tmp3 == "D" & tmp2$ratio <1) )
        tmp$keep[i] <- F
    }
    tmp <- tmp[tmp$keep,]
    cnv_by_gene <- cnv_by_gene[!cnv_by_gene$Segments %in% tmp$Segments,]
    
    write.csv(tmp[,-ncol(tmp)], "SCNA/important_CNVs.csv", row.names = F)
  }
}

### CGC + TS/OG Concordant
if(any(cnv_by_gene$Gene %in% CGC_df$Gene.Symbol)) {
  tmp <- cnv_by_gene[cnv_by_gene$Gene %in% CGC_df$Gene.Symbol,]
  tmp <- tmp[tmp$CN != 2,]
  tmp$keep <- T
  for(i in 1:nrow(tmp)) {
    tmp2 <- tmp[i,]
    tmp3 <- CGC_df$Role.in.Cancer[CGC_df$Gene.Symbol == tmp2$Gene]
    if( !( grepl("onc",tmp3) & tmp2$ratio >1) & !(grepl("TSG",tmp3) & tmp2$ratio <1) )
      tmp$keep[i] <- F
  }
  tmp <- tmp[tmp$keep,]
  
  write.csv(tmp[,-ncol(tmp)], "SCNA/CGC_CNVs.csv", row.names = F)
}

### Broad SCNAs
broad <- cnv[!is.na(cnv$cytoband),-1]
broad <- broad[,c(15,2,9,11)]
write.csv(broad, "SCNA/broad.csv", row.names = FALSE)

# LOH - exomeCNV ---------------------------------------------------------------------
dir.create("LOH")
# read in loh file
loh <- read.csv("../ExomeCNV/LOH_regions.csv", header = TRUE)
loh <- loh[loh$difference > 0.4,]
loh <- loh[,c("chr", "position.start", "position.end", "normal_b", "tumor_b", "difference")]

#find genes that are overlapped by segment and segments that overlap genes
loh_genes_overlap <- annotate_genes(loh)
loh_cytb_overlap <- annotate_cytb(loh, cytobands_df)

loh$length <- loh$position.end - loh$position.start + 1
loh$genes <- sapply(loh_genes_overlap$Segment_genes, function(x) ifelse(length(x) > 100, "100",
                                                                             paste(x, collapse = ", ")))
loh$cytoband <- sapply(loh_cytb_overlap$Segment_cytbs, function(x) shorten_cyt_anno(x, cytobands_df))

loh <- loh[, c("genes", "length", "cytoband", setdiff(colnames(loh), c("genes", "cytoband", "length")))]
colnames(loh) <- c("Genes", "Length", "Cytoband", "Chr", "Start", "End", "N_BAF", "T_BAF", "Absolute_Diff")
write.csv(loh, "LOH/LOH.csv", row.names = F)

loh_by_gene <- data.frame(Gene = names(loh_genes_overlap$Gene_segments), 
                          Segment = NA, Absolute_Diff = NA)
loh$id <- paste0(loh$Chr,":", loh$Start, "-", loh$End)
gene_segs <- loh_genes_overlap$Gene_segments

for(i in 1:length(gene_segs)) {
  genes_segments <- gene_segs[[i]]
  idx <- match(genes_segments, loh$id)
  
  loh_by_gene$Segment[i] <- loh$id[idx[1]]
  loh_by_gene$Absolute_Diff[i] <- loh$Absolute_Diff[idx[1]]
  if(length(idx) > 1) {
    for(j in 2:length(idx)) {
      loh_by_gene <- rbind(loh_by_gene, 
                           data.frame(Gene = loh_by_gene$Gene[i],
                                      Segment = loh$id[idx[j]],
                                      Absolute_Diff = loh$Absolute_Diff[idx[j]]))
    }
  }
}

### CGC genes
if(any(loh_by_gene$Gene %in% CGC_df$Gene.Symbol))
{
  stmp <- loh_by_gene[loh_by_gene$Gene %in% CGC_df$Gene.Symbol,]
  write.csv(tmp, "LOH/CGC_loh.csv", row.names = F)
}