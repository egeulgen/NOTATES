############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                Integration of Alterations                #
#                   Ege Ulgen, Dec 2017                    #
############################################################

dir.create("./NOTATES/")
setwd("./NOTATES/")

# dir for data sources (same as the directory where script is located)
initial_options <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", initial_options[grep("--file=", initial_options)])
script_dir <- dirname(script_name)

# Necessary resources -----------------------------------------------------
# CGC
CGC_df <- read.csv(paste0(script_dir, "/CGC_dec15_17.csv"), stringsAsFactors = F)

# Curated Databases
dna_repair_df <- read.csv(paste0(script_dir, "/curated_dbs/DNA_damage_repair_22feb16.csv"), stringsAsFactors = F)
irinotecan_df <- read.csv(paste0(script_dir, "/curated_dbs/irinotecan_response_genes_22feb16.csv"),stringsAsFactors = F)
tmz_df <- read.csv(paste0(script_dir, "/curated_dbs/temozolamide_resistance_genes_22feb16.csv"),stringsAsFactors = F)
KEGG_df <- read.csv(paste0(script_dir, "/curated_dbs/important_KEGG_pws_3jan17.csv"),stringsAsFactors = F)

# Curated Alterations
curated_SNV <- read.csv(paste0(script_dir, "/curated_alterations/glioma_important_SNV_May25_17.csv"), stringsAsFactors = F)
curated_CNA <- read.csv(paste0(script_dir, "/curated_alterations/glioma_important_CNA_May25_17.csv"), stringsAsFactors = F)

# Germline Mutations ------------------------------------------------------
dir.create("./Germline")
germline_mutations <- read.csv("../Germline/output/germline_variant_report.csv", stringsAsFactors = F)

### ACMG incidental
if(any(grepl("ACMG",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("ACMG",germline_mutations$Filter_Group),]
  tmp <- tmp[,c("Hugo_Symbol","id","Filter_Comment","Variant_Classification",
                "minor_allele","Reference_Allele", "Germline_Seq_Allele2", "allele_frequency")]
  colnames(tmp) <- c("Gene","ID","Disease(s)","Effect","Minor Allele", "Ref", "Var", "AF")
  germline_mutations <- germline_mutations[!grepl("ACMG",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/ACMG.csv", row.names = F)
}
### Cancer Gene Census
if(any(grepl("CGC",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("CGC",germline_mutations$Filter_Group),]
  tmp <- tmp[,c("Hugo_Symbol","id","Variant_Classification",
                "minor_allele","Reference_Allele", "Germline_Seq_Allele2", "allele_frequency")]
  colnames(tmp) <- c("Gene","ID","Effect","Minor Allele", "Ref", "Var", "AF")
  germline_mutations <- germline_mutations[!grepl("CGC",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/CGC.csv", row.names = F)
}
### Cancer Predisposition Gene
if(any(grepl("CPG",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("CPG",germline_mutations$Filter_Group),]
  tmp <- tmp[,c("Hugo_Symbol","id","Variant_Classification",
                "minor_allele","Reference_Allele", "Germline_Seq_Allele2", "allele_frequency")]
  colnames(tmp) <- c("Gene","ID","Effect","Minor Allele", "Ref", "Var", "AF")
  germline_mutations <- germline_mutations[!grepl("CPG",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/CPG.csv", row.names = F)
}
### Fanconi Anemia Pathway
if(any(grepl("FAP",germline_mutations$Filter_Group)))
{
  tmp <- germline_mutations[grepl("FAP",germline_mutations$Filter_Group),]
  tmp <- tmp[,c("Hugo_Symbol","id","Variant_Classification",
                "minor_allele","Reference_Allele", "Germline_Seq_Allele2", "allele_frequency")]
  colnames(tmp) <- c("Gene","ID","Effect","Minor Allele", "Ref", "Var", "AF")
  germline_mutations <- germline_mutations[!grepl("FAP",germline_mutations$Filter_Group),]
  
  write.csv(tmp, "./Germline/FAP.csv", row.names = F)
}
### Other
if(nrow(germline_mutations) != 0)
{
  tmp <- germline_mutations
  tmp <- tmp[,c("Hugo_Symbol","id","Variant_Classification",
                "minor_allele","Reference_Allele", "Germline_Seq_Allele2", "allele_frequency")]
  colnames(tmp) <- c("Gene","ID","Effect","Minor Allele", "Ref", "Var", "AF")
  
  write.csv(tmp, "./Germline/OTHER.csv", row.names = F)
}

### Common Variants
common_var <- read.csv("../Germline/output/common_variant_report.csv", stringsAsFactors = F)
common_var <- common_var[,c("id","Hugo_Symbol","Variant_Classification",
                            "gCCV_Risk_allele","Reference_Allele","Normal_Seq_Allele2","allele_frequency")]
colnames(common_var) <-  c("ID","Gene","Effect","Risk Allele", "Ref", "Var", "AF")
write.csv(common_var, "./Germline/common_var.csv", row.names = F)

# Somatic Mutations -------------------------------------------------------
dir.create("./Somatic_SNV")
somatic_SNVs <- read.delim("../Oncotator/annotated.sSNVs.tsv", stringsAsFactors=F, comment.char="#")

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
    tmp <- somatic_SNVs[i,]
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
  write.csv(tmp, "Somatic_SNV/important_glioma_SNVs.csv",row.names = F)
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

# Helper functions for annotation and Cytoband & Gene annotations ---------
find_feat_in_seg <- function(seg_df, feat_df, 
                             chrom1 = 1, start1 = 2, end1 = 3, 
                             chrom2 = 1, start2 = 2, end2 = 3, feature = 4, threshold = .25) {
  # find overlaps
  chr.list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  feat_segments <- list()
  Seg_feats <- list()
  
  feat_df$width <- feat_df[, end2] - feat_df[, start2] + 1
  
  for(chrom in chr.list)
  {
    tmpc <- seg_df[seg_df[, chrom1] == chrom, ]
    tempf <- feat_df[feat_df[, chrom2] == chrom, ]
    
    if(nrow(tmpc) != 0)
    {
      for(i in 1:nrow(tmpc))
      {
        tmp_feats <- c()
        Ss <- tmpc[i, start1]
        Se <- tmpc[i, end1]
        
        idx <- which(tempf[, start2] < Se & tempf[, end2] > Ss)
        if(length(idx) != 0)
        {
          for(j in idx)
          {
            Fs <- tempf[j, start2]
            Fe <- tempf[j, end2]
            # calculate the proportion
            if(Fs >= Ss & Fe <= Se) # segment contains gene
              prop <- 1
            else if(Ss >= Fs & Se <= Fe) # gene contains segment
              prop <- (Se - Ss + 1)/tempf$width[j]
            else if(Fs < Ss & Fe < Se) # gene intersects segment on left
              prop <- (Fe - Ss + 1)/tempf$width[j]
            else if(Fe > Se & Fs > Ss) # gene intersects segment on right
              prop <- (Se - Fs + 1)/tempf$width[j]
            
            ## if prop >=threshold(default is .25), accept that the segment contains the gene
            if( prop >= threshold )
              tmp_feats <- c(tmp_feats, tempf[j, feature])
            
            ## add to features list
            names(prop) <- paste0(chrom, ":", Ss, "-", Se)
            feat_segments[[tempf[j, feature]]] <- c(feat_segments[[tempf[j, feature]]], prop)
          }
        }
        ## add segment's genes to the list
        Seg_feats <- append(Seg_feats, list(tmp_feats))
      }
    }
  }
  return(list(Feat_segments = feat_segments, Segment_feats = Seg_feats))
}

cyt_func <- function(cyt, cytobands_df)
{
  if(length(cyt) == 1) ## single cytoband
    return(cyt)
  else if(is.null(cyt))
    return(NA)
  else 
  {
    tmp <- subset(cytobands, V1 == sub("[p-q].*", "", cyt[1]))
    
    tmp_P <- tmp$V4[grepl("p", tmp$V4)]
    tmp_P_start <- tmp_P[1]
    tmp_P_end <- tmp_P[length(tmp_P)]
    
    tmp_Q <- tmp$V4[grepl("q", tmp$V4)]
    tmp_Q_start <- tmp_Q[1]
    tmp_Q_end <- tmp_Q[length(tmp_Q)]
    
    if (cyt[1] == tmp_P_start & cyt[length(cyt)] == tmp_Q_end) ## whole chromosome
      return(sub("[p-q].*", "", cyt[1]))
    else if (cyt[1] == tmp_P_start & cyt[length(cyt)] == tmp_P_end) ## whole short arm
      return(paste0(sub("[p-q].*", "", cyt[1]), "p"))
    else if (cyt[1] == tmp_Q_start & cyt[length(cyt)] == tmp_Q_end) ## whole long arm
      return(paste0(sub("[p-q].*", "", cyt[1]), "q"))
    else ## multiple segments
      return(paste0(cyt[1], "-", cyt[length(cyt)]))
  }
}

# read in genes file
HS_genes <- read.csv(paste0(script_dir, "/homo_sapiens_genes.csv"), stringsAsFactors = F)
HS_genes <- HS_genes[,-c(4,5)]

# read in cytobands file
cytobands <- read.delim(paste0(script_dir, "/cytoBand.txt"), header = F, stringsAsFactors = F)
cytobands$V4 <- paste0(cytobands$V1, cytobands$V4)

# CNV - exomeCNV ---------------------------------------------------------------------
dir.create("./SCNA")
# read in cnv file
cnv <- read.delim("../ExomeCNV/CNV.cnv.txt", header = T, stringsAsFactors = F)
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

cnv$copy.number <- sapply(cnv$ratio, copy_num_call)

#find genes that are overlapped by segment and segments that overlap genes
cnv_genes_overlap <- find_feat_in_seg(cnv, HS_genes)
cnv_cytb_overlap <- find_feat_in_seg(cnv, cytobands, threshold = .5)

cnv$length <- cnv$probe_end - cnv$probe_start + 1
cnv$genes <- sapply(cnv_genes_overlap[["Segment_feats"]], function(x) ifelse(length(x) > 500, ">500", paste(x, collapse = ", ")))
cnv$cytoband <- sapply(cnv_cytb_overlap[["Segment_feats"]], function(cyt) cyt_func(cyt, cytobands))

cnv <- cnv[, c("genes", "length", "cytoband", setdiff(colnames(cnv), c("genes", "cytoband", "length")))]
write.csv(cnv, "SCNA/exomeCNV.csv", row.names = F)

cnv_by_gene <- data.frame(Gene=names(cnv_genes_overlap[["Feat_segments"]]),
                          Segment = NA, 
                          ratio = NA, CN = NA, av_cov = NA, stringsAsFactors = F)
cnv$id <- paste0(cnv$chr,":", cnv$probe_start, "-", cnv$probe_end)
gene_segs <- cnv_genes_overlap[["Feat_segments"]]
N <- length(gene_segs)
for(i in 1:N)
{
  gene <- gene_segs[[i]]
  idx <- match(names(gene), cnv$id)
  
  tmp_ratio <- cnv$ratio[idx]
  tmp_av_cov <- cnv$average.coverage[idx]
  
  if(length(gene) == 1)
  {
    cnv_by_gene$ratio[i] <- cnv$ratio[idx]
    cnv_by_gene$av_cov[i] <- cnv$average.coverage[idx]
    cnv_by_gene$Segment[i] <- names(gene)
  }else
  {
    cnv_by_gene$ratio[i] <- cnv$ratio[idx[1]]
    cnv_by_gene$av_cov[i] <- cnv$average.coverage[idx[1]]
    cnv_by_gene$Segment[i] <- names(gene)[1]
    for(j in 2:length(gene))
    {
      cnv_by_gene <- rbind(cnv_by_gene, c(names(gene_segs)[i], 
                                          names(gene)[j],
                                          cnv$ratio[idx[j]],
                                          NA,
                                          cnv$average.coverage[idx[j]]))
    }
  }
}
cnv_by_gene$av_cov <- as.numeric(cnv_by_gene$av_cov)
cnv_by_gene$ratio <- as.numeric(cnv_by_gene$ratio)
cnv_by_gene$CN <- sapply(cnv_by_gene$ratio, copy_num_call)

write.csv(cnv_by_gene, "SCNA/all_genes.csv", row.names = F)

### Glioma important SCNAs
if(any(cnv_by_gene$Gene %in% curated_CNA$Gene))
{
  tmp <- cnv_by_gene[cnv_by_gene$Gene %in% curated_CNA$Gene,]
  tmp <- tmp[tmp$CN != 2,]
  tmp$keep <- T
  for(i in 1:nrow(tmp))
  {
    tmp2 <- tmp[i,]
    tmp3 <- curated_CNA$Type[curated_CNA$Gene == tmp2$Gene]
    if( !(tmp3 == "A" & tmp2$ratio >1) & !(tmp3 == "D" & tmp2$ratio <1) )
      tmp$keep[i] <- F
  }
  tmp <- tmp[tmp$keep,]
  cnv_by_gene <- cnv_by_gene[!cnv_by_gene$Segments %in% tmp$Segments,]
    
  write.csv(tmp[,-ncol(tmp)], "SCNA/important_CNVs.csv", row.names = F)
}

### CGC + TS/OG Concordant
if(any(cnv_by_gene$Gene %in% CGC_df$Gene.Symbol))
{
  tmp <- cnv_by_gene[cnv_by_gene$Gene %in% CGC_df$Gene.Symbol,]
  tmp <- tmp[tmp$CN != 2,]
  tmp$keep <- T
  for(i in 1:nrow(tmp))
  {
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
write.csv(broad, "SCNA/broad.csv", row.names = F)

# LOH - exomeCNV ---------------------------------------------------------------------
dir.create("LOH")
# read in loh file
loh <- read.csv("../ExomeCNV/LOH_regions.csv", header = T, stringsAsFactors = F)
loh$LOH <- NULL

#find genes that are overlapped by segment and segments that overlap genes
loh_genes_overlap <- find_feat_in_seg(loh, HS_genes)
loh_cytb_overlap <- find_feat_in_seg(loh, cytobands, threshold = .5)

loh$length <- loh$position.end - loh$position.start + 1
loh$genes <- sapply(loh_genes_overlap[["Segment_feats"]], function(x) ifelse(length(x) > 500, ">500", paste(x, collapse = ", ")))
loh$cytoband <- sapply(loh_cytb_overlap[["Segment_feats"]], function(cyt) cyt_func(cyt, cytobands))

loh <- loh[, c("genes", "length", "cytoband", setdiff(colnames(loh), c("genes", "cytoband", "length")))]
write.csv(loh, "LOH/LOH.csv", row.names = F)

loh_by_gene <- data.frame(Gene=names(loh_genes_overlap[["Feat_segments"]]), 
                          num_segments = NA, LOH = T, stringsAsFactors = F)
loh$id <- paste0(loh$chr,":", loh$position.start, "-", loh$position.end)
gene_segs <- loh_genes_overlap[["Feat_segments"]]
for(i in 1:length(gene_segs))
{
  gene <- gene_segs[[i]]
  idx <- match(names(gene), loh$id)
  loh_by_gene$num_segments[i] <- length(gene)
}

### CGC genes
if(any(loh_by_gene$Gene %in% CGC_df$Gene.Symbol))
{
  tmp <- loh_by_gene[loh_by_gene$Gene %in% CGC_df$Gene.Symbol,]
  write.csv(tmp, "LOH/CGC_loh.csv", row.names = F)
}
