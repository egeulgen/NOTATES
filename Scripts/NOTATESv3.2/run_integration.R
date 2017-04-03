############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                Integration of Alterations                #
#                   Ege Ulgen, Mar 2017                    #
############################################################

dir.create("./NOTATES/")
setwd("./NOTATES/")

# dir for data sources (same as the directory where script is located)
initial_options <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", initial_options[grep("--file=", initial_options)])
script_dir <- dirname(script_name)

# Necessary resources -----------------------------------------------------
# TCGA tables for alteration percentages
tcga_idh_MUT <- paste0(script_dir, "/TCGA_feb21_16/TCGA_IDH_MUT_df.csv")
tcga_idh_WT <- paste0(script_dir, "/TCGA_feb21_16/TCGA_IDH_WT_df.csv")

# CGC
CGC_df <- read.csv(paste0(script_dir, "/cancer_gene_census.csv"), stringsAsFactors = F)

# KEGG Patways in Cancer genes
kegg_pathways_in_cancer <- read.delim(paste0(script_dir, "/KEGG_pathways_in_cancer.txt"), stringsAsFactors = F)
kegg_pathways_in_cancer <- kegg_pathways_in_cancer[,1]

# Curated Databases
dna_repair_df <- read.csv(paste0(script_dir, "/curated_dbs/DNA_damage_repair_22feb16.csv"), stringsAsFactors = F)
irinotecan_df <- read.csv(paste0(script_dir, "/curated_dbs/irinotecan_response_genes_22feb16.csv"),stringsAsFactors = F)
tmz_df <- read.csv(paste0(script_dir, "/curated_dbs/temozolamide_resistance_genes_22feb16.csv"),stringsAsFactors = F)
KEGG_df <- read.csv(paste0(script_dir, "/curated_dbs/important_KEGG_pws_3jan17.csv"),stringsAsFactors = F)

# Somatic Mutations -------------------------------------------------------
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

# Reject certain non-coding variant types not likely to be interpreted as important
somatic_SNVs <- subset(somatic_SNVs,
                       !( Variant_Classification %in% c("Silent", "3'UTR", "3'Flank", "5'UTR", 
                                                        "5'Flank", "IGR", "Intron", "lincRNA", "RNA") ))

# Reject variants with "Unknown" Hugo Symbols
somatic_SNVs <- subset(somatic_SNVs, Hugo_Symbol != "Unknown")

# Change "" protein changes to NA
somatic_SNVs$Protein_Change[somatic_SNVs$Protein_Change == ""] <- NA

# keep only needed annotations
keep <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", "Variant_Type", 
          "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "tumor_f", "genotype", "dbSNP_RS", "Genome_Change", "Codon_Change",
          "Protein_Change", "UniProt_AApos", "DNARepairGenes_Role", "FamilialCancerDatabase_Syndromes",
          "COSMIC_n_overlapping_mutations", "COSMIC_total_alterations_in_gene", 
          colnames(somatic_SNVs)[grepl("^GO_", colnames(somatic_SNVs))], 
          colnames(somatic_SNVs)[grepl("^CGC_", colnames(somatic_SNVs))])

write.csv(somatic_SNVs[, keep], "detailed_sSNV.csv", row.names = F, quote = F)

somatic_SNVs <- somatic_SNVs[, c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "Genome_Change",
                                 "tumor_f", "genotype", "COSMIC_n_overlapping_mutations",
                                 "COSMIC_total_alterations_in_gene")]

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
# read in cnv file
cnv <- read.delim("../ExomeCNV/CNV.cnv.txt", header = T, stringsAsFactors = F)
# discard rows with NA copy.number
cnv <- subset(cnv, !is.na(copy.number))
# discard rows with NA sens
cnv <- subset(cnv, !is.na(sens))

# mean and sem for cutoff
cutoff <- c()
cutoff$mean <- mean(cnv$ratio)
cutoff$sem <- sd(cnv$ratio)/sqrt(length(cnv$ratio))

copy_num_call <- function(ratio, cutoff){
  if( ( ratio >= (cutoff$mean -3*cutoff$sem) ) & ( ratio <= (cutoff$mean + 3*cutoff$sem) ) )
    return(2)
  else if( ratio < cutoff$mean - 3*cutoff$sem & ratio >= 0.5 )
    return(1)
  else if( ratio < 0.5 )
    return(0)
  else
    return(ceiling(ratio*2))
}

cnv$copy.number <- sapply(cnv$ratio, function(x) copy_num_call(x, cutoff))

#find genes that are overlapped by segment and segments that overlap genes
cnv_genes_overlap <- find_feat_in_seg(cnv, HS_genes)
cnv_cytb_overlap <- find_feat_in_seg(cnv, cytobands, threshold = .5)

cnv$length <- cnv$probe_end - cnv$probe_start + 1
cnv$genes <- sapply(cnv_genes_overlap[["Segment_feats"]], function(x) ifelse(length(x) > 500, ">500", paste(x, collapse = ", ")))
cnv$cytoband <- sapply(cnv_cytb_overlap[["Segment_feats"]], function(cyt) cyt_func(cyt, cytobands))

cnv <- cnv[, c("genes", "length", "cytoband", setdiff(colnames(cnv), c("genes", "cytoband", "length")))]
write.csv(cnv, "exomeCNV.csv", row.names = F)


cnv_by_gene <- data.frame(Gene=names(cnv_genes_overlap[["Feat_segments"]]), 
                          num_segments = NA,
                          Segments = NA, 
                          ratio = NA, CN = NA, av_cov = NA, stringsAsFactors = F)
cnv$id <- paste0(cnv$chr,":", cnv$probe_start, "-", cnv$probe_end)
gene_segs <- cnv_genes_overlap[["Feat_segments"]]
for(i in 1:length(gene_segs))
{
  gene <- gene_segs[[i]]
  idx <- match(names(gene), cnv$id)
  tmp_ratio <- cnv$ratio[idx]
  tmp_av_cov <- cnv$average.coverage[idx]
  
  cnv_by_gene$num_segments[i] <- length(gene)
  
  cnv_by_gene$ratio[i] <- sum(tmp_ratio * gene)/sum(gene)
  cnv_by_gene$av_cov[i] <- paste(round(tmp_av_cov,2), collapse = ", ") 
  cnv_by_gene$Segments[i] <- paste(names(gene), collapse = ", ")
}
cnv_by_gene$CN <- sapply(cnv_by_gene$ratio, function(x) copy_num_call(x, cutoff))




# LOH - exomeCNV ---------------------------------------------------------------------
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
write.csv(loh, "LOH.csv", row.names = F)

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


# CNV - THetA -----------------------------------------------------------
theta_segments <- read.delim("../THetA/output/CNV.n2.withBounds", stringsAsFactors = F)
theta_result <- read.delim("../THetA/output/CNV.n2.results", header = T, stringsAsFactors = F)

theta_segments$CopyNumber <- unlist(strsplit(theta_result$C, ":"))
theta_segments$CopyNumber[theta_segments$CopyNumber == "X"] <- NA
theta_segments$CopyNumber <- as.numeric(theta_segments$CopyNumber)

theta_segments$chrm <- paste0("chr", theta_segments$chrm)
theta_segments$chrm <- gsub("chr23", "chrX", theta_segments$chrm)
theta_segments$chrm <- gsub("chr24", "chrY", theta_segments$chrm)

theta_segments <- theta_segments[!is.na(theta_segments$CopyNumber), ]
theta_segments <- theta_segments[,-1]

#find genes that are overlapped by segment and segments that overlap genes
theta_genes_overlap <- find_feat_in_seg(theta_segments, HS_genes)
theta_cytb_overlap <- find_feat_in_seg(theta_segments, cytobands, threshold = .5)

theta_segments$length <- theta_segments$end - theta_segments$start + 1
theta_segments$genes <- sapply(theta_genes_overlap[["Segment_feats"]], function(x) ifelse(length(x) > 500, ">500", paste(x, collapse = ", ")))
theta_segments$cytoband <- sapply(theta_cytb_overlap[["Segment_feats"]], function(cyt) cyt_func(cyt, cytobands))

theta_segments <- theta_segments[, c("genes", "length", "cytoband", setdiff(colnames(theta_segments), c("genes", "cytoband", "length")))]
write.csv(theta_segments, "THetA.csv", row.names = F)

theta_by_gene <- data.frame(Gene=names(theta_genes_overlap[["Feat_segments"]]), 
                            num_segments = NA,
                            Segments=NA, 
                            CN=NA, stringsAsFactors = F)
theta_segments$id <- paste0(theta_segments$chr,":", theta_segments$start, "-", theta_segments$end)
gene_segs <- theta_genes_overlap[["Feat_segments"]]
for(i in 1:length(gene_segs))
{
  gene <- gene_segs[[i]]
  idx <- match(names(gene), theta_segments$id)
  tmp_cn <- theta_segments$CopyNumber[idx]
  
  theta_by_gene$num_segments[i] <- length(gene)
  theta_by_gene$CN[i] <- round(sum(tmp_cn * gene)/sum(gene))
  theta_by_gene$Segments[i] <- paste(names(gene), collapse = ", ")
}

# Germline Mutations ------------------------------------------------------
germline_mutations <- read.csv("../Germline/output/germline_variant_report.csv", stringsAsFactors = F)

# Final Table -------------------------------------------------------------
#create the table
genes <- unique(c(somatic_SNVs$Hugo_Symbol, germline_mutations$Hugo_Symbol))
alteration_table <- as.data.frame(genes, stringsAsFactors = F)
colnames(alteration_table) <- "Hugo_Symbol"

#choose the appropriate TCGA data between IDH-mutant or IDH-wildtype group
if(any(somatic_SNVs$Hugo_Symbol %in% c("IDH1","IDH2"))){
  tcga_df <- read.csv(tcga_idh_MUT, stringsAsFactors = F)
} else {
  tcga_df <- read.csv(tcga_idh_WT, stringsAsFactors = F) 
}

#### Somatic SNVs
alteration_table$sSNV_number <- 0
alteration_table$sSNV_class <- NA
alteration_table$sSNV_Genome_Change <- NA
alteration_table$sSNV_Protein_Change <- NA
alteration_table$sSNV_allele_frequency <- NA
alteration_table$sSNV_n_overlapping <- NA
alteration_table$sSNV_total_COSMIC<- NA
for(i in which(alteration_table$Hugo_Symbol %in% somatic_SNVs$Hugo_Symbol))
{
  tmp <- somatic_SNVs[somatic_SNVs$Hugo_Symbol == alteration_table$Hugo_Symbol[i], ]
  
  alteration_table$sSNV_number[i] <- nrow(tmp)
  alteration_table$sSNV_class[i] <- paste(tmp$Variant_Classification, collapse = ", ")
  alteration_table$sSNV_Genome_Change[i] <- paste(tmp$Genome_Change, collapse = ", ")
  alteration_table$sSNV_Protein_Change[i] <- paste(tmp$Protein_Change, collapse = ", ")
  alteration_table$sSNV_allele_frequency[i] <- paste(tmp$tumor_f, collapse = ", ")
  alteration_table$sSNV_n_overlapping[i] <- paste(tmp$COSMIC_n_overlapping_mutations, collapse = ", ")
  alteration_table$sSNV_total_COSMIC[i] <- paste(tmp$COSMIC_total_alterations_in_gene, collapse = ", ")
}

#### Germline SNPs
alteration_table$Germline_SNP <- 0
alteration_table$Germline_RS_ID <- NA
alteration_table$Germline_Filter_Group <- NA
alteration_table$Germline_Frequency <- NA
for(i in which(alteration_table$Hugo_Symbol %in% germline_mutations$Hugo_Symbol))
{
  tmp <- germline_mutations[germline_mutations$Hugo_Symbol == alteration_table$Hugo_Symbol[i], ]
  alteration_table$Germline_SNP[i] <- nrow(tmp)
  alteration_table$Germline_RS_ID[i] <- paste(tmp$id, collapse = ", ")
  alteration_table$Germline_Filter_Group[i] <- paste(tmp$Filter_Group, collapse = ", ")
  alteration_table$Germline_Frequency[i] <- paste(tmp$allele_frequency, collapse = ", ")
}

#### CNV - ExomeCNV
## add CNA info to genes with single-nucleotide variations
alteration_table$ExomeCNV_CN <- NA
alteration_table$ExomeCNV_ratio <- NA
alteration_table$ExomeCNV_ave_coverage <- NA
alteration_table$ExomeCNV_num_seg <- NA
alteration_table$ExomeCNV_chr_pos <- NA

if(any(alteration_table$Hugo_Symbol %in% cnv_by_gene$Gene))
{
  for(i in which(alteration_table$Hugo_Symbol %in% cnv_by_gene$Gene))
  {
    idx <- which(cnv_by_gene$Gene == alteration_table$Hugo_Symbol[i])
    alteration_table$ExomeCNV_CN[i] <- cnv_by_gene$CN[idx]
    alteration_table$ExomeCNV_ratio[i] <- cnv_by_gene$ratio[idx]
    alteration_table$ExomeCNV_ave_coverage <- cnv_by_gene$av_cov[idx]
    alteration_table$ExomeCNV_chr_pos[i] <- cnv_by_gene$Segments[idx]
    alteration_table$ExomeCNV_num_seg[i] <- cnv_by_gene$num_segments[idx]
  }
}

## add additional genes with CNA only, only KEGG Pathways in Cancer genes
idx <- which(cnv_by_gene$Gene %in% kegg_pathways_in_cancer)
template <- rep(NA, ncol(alteration_table))
names(template) <- colnames(alteration_table)

for(i in idx)
{
  gene <- cnv_by_gene$Gene[i]
  if(!(gene %in% alteration_table$Hugo_Symbol))
  {
    tmp <- template
    tmp["Hugo_Symbol"] <- gene
    tmp["sSNV_number"] <- 0
    tmp["Germline_SNP"] <- 0
    tmp["ExomeCNV_CN"] <- cnv_by_gene$CN[i]
    tmp["ExomeCNV_ratio"] <- cnv_by_gene$ratio[i]
    tmp["ExomeCNV_ave_coverage"] <- cnv_by_gene$av_cov[i]
    tmp["ExomeCNV_num_seg"] <- cnv_by_gene$num_segments[i]
    tmp["ExomeCNV_chr_pos"] <- cnv_by_gene$Segments[i]
    
    alteration_table <- rbind(alteration_table, tmp)
  }
}

## CNA segments
cnv$cytoband[is.na(cnv$cytoband)] <- cnv$id[is.na(cnv$cytoband)]

idx <- c()
for(i in 1:length(cnv$id))
{
  if(!any(grepl(cnv$id[i], alteration_table$ExomeCNV_chr_pos)))
    idx <- c(idx, i)
}

for(i in idx)
{
  tmp <- template
  tmp["ExomeCNV_CN"] <- cnv$copy.number[i]
  tmp["ExomeCNV_ratio"] <- cnv$ratio[i]
  tmp["ExomeCNV_ave_coverage"] <- cnv$average.coverage[i]
  tmp["ExomeCNV_num_seg"] <- 1
  tmp["Hugo_Symbol"] <- cnv$cytoband[i]
  tmp["ExomeCNV_chr_pos"] <- cnv$id[i]
  
  alteration_table <- rbind(alteration_table, tmp)
}

#### CNV - THetA
# genes
idx <- which(alteration_table$Hugo_Symbol %in% theta_by_gene$Gene)
alteration_table$THetaCN <- NA
for(i in idx)
{
  gene <- alteration_table$Hugo_Symbol[i]
  alteration_table$THetaCN[i] <- theta_by_gene$CN[theta_by_gene$Gene == gene]
}
# regions
idx <- which(alteration_table$Hugo_Symbol %in% theta_segments$id)
for(i in idx)
{
  region <- alteration_table$Hugo_Symbol[i]
  alteration_table$THetaCN[i] <- theta_segments$CopyNumber[theta_segments$id == region]
}

#### LOH - ExomeCNV
# genes
idx <- which(alteration_table$Hugo_Symbol %in% loh_by_gene$Gene)
alteration_table$LOH <- NA
for(i in idx)
{
  gene <- alteration_table$Hugo_Symbol[i]
  alteration_table$LOH[i] <- loh_by_gene$LOH[loh_by_gene$Gene == gene]
}

## add additional genes with LOH only, only KEGG Pathways in Cancer genes
idx <- which(loh_by_gene$Gene %in% kegg_pathways_in_cancer)
template <- rep(NA, ncol(alteration_table))
names(template) <- colnames(alteration_table)

for(i in idx)
{
  gene <- loh_by_gene$Gene[i]
  if(!(gene %in% alteration_table$Hugo_Symbol))
  {
    tmp <- template
    tmp["Hugo_Symbol"] <- gene
    tmp["sSNV_number"] <- 0
    tmp["Germline_SNP"] <- 0
    tmp["LOH"] <- TRUE
    alteration_table <- rbind(alteration_table, tmp)
  }
}

#### CGC
alteration_table$inCGC <- F
alteration_table$CGC_Genetics <- NA
alteration_table$CNA_concordant_w_CGC <- F
for(i in which(alteration_table$Hugo_Symbol %in% CGC_df$Gene.Symbol))
{
  alteration_table$inCGC[i] <- T
  tmp <- CGC_df[CGC_df$Gene.Symbol == alteration_table$Hugo_Symbol[i], ]
  
  # CGC molecular genetics
  alteration_table$CGC_Genetics[i] <- tmp$Molecular.Genetics
  
  # CNA concordant with CGC molecular genetics?
  tmp <- alteration_table$THetaCN[i]
  cond1 <- (tmp < 2 && grepl("Rec", alteration_table$CGC_Genetics[i]))
  cond2 <- (tmp > 2 && grepl("Dom", alteration_table$CGC_Genetics[i]))
  if(cond1 || cond2)
    alteration_table$CNA_concordant_w_CGC[i] = T
}
#### Clinical
alteration_table$Clinical <- alteration_table$inCGC & alteration_table$CNA_concordant_w_CGC
tmp <- strsplit(alteration_table$sSNV_n_overlapping,", ",fixed = T)
tmp <- sapply(tmp, as.numeric)
tmp <- sapply(tmp, function(x) any(na.omit(x) > 0))
alteration_table$Clinical <- (alteration_table$inCGC & tmp) | alteration_table$Clinical

#### DNA Damage Repair List
alteration_table$DNA_damage_repair <- NA
for(i in which(alteration_table$Hugo_Symbol %in% dna_repair_df$Gene.Name))
{
  tmp <- dna_repair_df[dna_repair_df$Gene.Name == alteration_table$Hugo_Symbol[i], ]
  alteration_table$DNA_damage_repair[i] <- tmp$FUNCTION
}

#### TMZ response
alteration_table$TMZ_response_corr <- NA
for(i in which(alteration_table$Hugo_Symbol %in% tmz_df$gene))
{
  tmp <- tmz_df[tmz_df$gene == alteration_table$Hugo_Symbol[i], ]
  alteration_table$TMZ_response_corr[i] <- tmp$correlation.w.temozolamide.response
}

#### Irinotecan response
alteration_table$Irinotecan_response <- alteration_table$Hugo_Symbol %in% irinotecan_df$Genes

#### KEGG pathways
alteration_table$Cell_Cycle <- F
alteration_table$mTOR_Signalling <- F
alteration_table$Pathways_in_Cancer <- F
for(i in which(alteration_table$Hugo_Symbol %in% KEGG_df$Gene))
{
  tmp <- KEGG_df[KEGG_df$Gene == alteration_table$Hugo_Symbol[i], ]
  
  if(any(tmp$pathway == "Cell Cycle"))
    alteration_table$Cell_Cycle[i] <- T
  if(any(tmp$pathway == "mTOR Signalling"))
    alteration_table$mTOR_Signalling[i] <- T
  if(any(tmp$pathway == "Pathways in Cancer"))
    alteration_table$Pathways_in_Cancer[i] <- T
}

#### TCGA
alteration_table$TCGA_mut_perc <- NA
alteration_table$TCGA_AMP_perc <- NA
alteration_table$TCGA_DEL_perc <- NA
alteration_table$TCGA_any_perc <- NA
for(i in which(alteration_table$Hugo_Symbol %in% tcga_df$Hugo_Symbol))
{
  tmp <- tcga_df[tcga_df$Hugo_Symbol == alteration_table$Hugo_Symbol[i], ]
  alteration_table$TCGA_mut_perc[i] <- tmp$Percent_Mut
  alteration_table$TCGA_AMP_perc[i] <- tmp$Percent_AMP
  alteration_table$TCGA_DEL_perc[i] <- tmp$Percent_DEL
  alteration_table$TCGA_any_perc[i] <- tmp$Percent_AnyAlteration
}

write.csv(alteration_table, "Alteration_table.csv", row.names = F)