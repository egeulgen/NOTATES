suppressPackageStartupMessages(library(StructuralVariantAnnotation))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(circlize))

# Read, filter and create necessary data ----------------------------------
# read VCF
vcf <- readVcf("DELLY/somatic_SV.vcf", "hg38")

# Extracting the structural variants as a GRanges
bpgr <- breakpointRanges(vcf, inferMissingBreakends = TRUE)

# filter for FILTER=="PASS" and svtype=="BND"
bpgr_filtered <- bpgr[bpgr$FILTER == "PASS" & bpgr$svtype == "BND", ]

# Convert breakpoint GRanges object to a Pairs object
pairs <- breakpointgr2pairs(bpgr_filtered)

# Filter for gene-to-different-gene SVs -----------------------------------
tra_df <- as.data.frame(pairs)

# Possible locations are ‘coding’, ‘intron’, ‘threeUTR’, ‘fiveUTR’, ‘intergenic’, ‘spliceSite’ and ‘promoter’.
valid_locs <- c("coding", "intron", "spliceSite")

if (nrow(tra_df) == 0) {
  tra_df <- cbind(tra_df, read.csv(text="Gene1,Gene2"))
} else {
  tra_df$Gene1 <- tra_df$Gene2 <- NA
  for (i in 1:nrow(tra_df)) {
    ### annotate first
    tmp_df_1 <- data.frame(chrom = tra_df$first.X.seqnames[i],
                           start = tra_df$first.X.start[i],
                           end = tra_df$first.X.end[i])
    locs_1 <- suppressWarnings(locateVariants(GenomicRanges::makeGRangesFromDataFrame(tmp_df_1), 
                                              TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                              AllVariants()))
    locs_1 <- locs_1[locs_1$LOCATION %in% valid_locs, ]
    locs_1 <- locs_1[!is.na(locs_1$GENEID), ]
    
    symbols_1 <- c()
    if (length(locs_1) != 0) {
      symbols_1 <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(locs_1$GENEID), "SYMBOL", "ENTREZID")
    }
    
    ### annotate second
    tmp_df_2 <- data.frame(chrom = tra_df$second.X.seqnames[i],
                           start = tra_df$second.X.start[i],
                           end = tra_df$second.X.end[i])
    locs_2 <- suppressWarnings(locateVariants(GenomicRanges::makeGRangesFromDataFrame(tmp_df_2), 
                                              TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                              AllVariants()))
    locs_2 <- locs_2[locs_2$LOCATION %in% valid_locs, ]
    locs_2 <- locs_2[!is.na(locs_2$GENEID), ]
    
    symbols_2 <- c()
    if (length(locs_2) != 0) {
      symbols_2 <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(locs_2$GENEID), "SYMBOL", "ENTREZID")
    }
    
    if (length(symbols_1) != 0 & length(symbols_2) != 0) {
      
      symbols_1 <- paste(symbols_1, collapse = "; ")
      symbols_2 <- paste(symbols_2, collapse = "; ")
      
      tra_df$Gene1[i] <- symbols_1
      tra_df$Gene2[i] <- symbols_2
      
    }
  }
}

# Visualize as Circos -----------------------------------------------------
final_pairs <- pairs[!is.na(tra_df$Gene1) & !is.na(tra_df$Gene2), ]
tra_final_df <- tra_df[!is.na(tra_df$Gene1) & !is.na(tra_df$Gene2), ]

### for gene labels
tmp_df <- tra_final_df[, c("first.X.seqnames", "first.X.start", "first.X.end", "Gene1")]
colnames(tmp_df) <- c("second.X.seqnames", "second.X.start", "second.X.end", "Gene2")
tmp_df <- rbind(tmp_df,
                tra_final_df[, c("second.X.seqnames", "second.X.start", "second.X.end", "Gene2")])

pdf("DELLY/TR_circos.pdf", width = 6, height = 6)
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(as.data.frame(S4Vectors::first(final_pairs)), 
                   as.data.frame(S4Vectors::second(final_pairs)))
if (nrow(tmp_df) != 0)
  circos.genomicLabels(tmp_df, labels.column = 4, side = "inside")
dev.off()

# create table of possible translocations ---------------------------------
final_table <- tra_final_df[, c("Gene1", "first.X.seqnames", "first.X.start", "first.X.end",
                                "Gene2", "second.X.seqnames", "second.X.start", "second.X.end")]
colnames(final_table) <- c("Gene1", "Chrom1", "Start1", "End1",
                           "Gene2", "Chrom2", "Start2", "End2")

write.csv(final_table, "DELLY/possible_translocations.csv", row.names = FALSE)

