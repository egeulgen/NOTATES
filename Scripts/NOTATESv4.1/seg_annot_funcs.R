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
            
            ## if prop >= threshold(default is .25), accept that the segment contains the gene
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