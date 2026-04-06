#!/usr/bin/env Rscript

library(dplyr)
library(circlize)

#######################################
# Logging setup
#######################################

log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("Starting Circos Plot generation\n")

#######################################
# Input files
#######################################

arriba_result <- snakemake@input[["arriba_result"]]


#######################################
# Output files
#######################################

circos_plot <- snakemake@output[["circos_plot"]]

#######################################
# Read results
#######################################
arriba_df <- read.delim(arriba_result)

#######################################
# Process
#######################################

arriba_df <- arriba_df %>%
  mutate(
    chr1 = paste0("chr", sub(":.*", "", breakpoint1)),
    pos1 = as.numeric(sub(".*:", "", breakpoint1)),
    chr2 = paste0("chr", sub(":.*", "", breakpoint2)),
    pos2 = as.numeric(sub(".*:", "", breakpoint2))
  )
arriba_plot <- arriba_df %>% filter(confidence == "high")

type_colors <- c(
  duplication = "#d7191c",
  translocation = "#fdae61",
  inversion = "#abdda4",
  deletion = "#2b83ba"
)

arriba_plot$color <- ifelse(
  grepl("duplication", arriba_plot$type),
  type_colors["duplication"],
  ifelse(
    grepl("translocation", arriba_plot$type),
    type_colors["translocation"],
    ifelse(
      grepl("inversion", arriba_plot$type),
      type_colors["inversion"],
      type_colors["deletion"]
    )
  )
)

#######################################
# Create Circos Plot
#######################################

cat("Creating Circos Plot...\n")

pdf(circos_plot, width = 8, height = 8)

circos.clear()
circos.par("start.degree" = 90, gap.degree = 2)
# Initialize chromosomes (human hg38)
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
circos.initializeWithIdeogram(species = "hg38", chromosome.index = chromosomes)


for (i in 1:nrow(arriba_plot)) {
  circos.link(
    sector.index1 = arriba_plot$chr1[i],
    point1 = arriba_plot$pos1[i],
    sector.index2 = arriba_plot$chr2[i],
    point2 = arriba_plot$pos2[i],
    col = ifelse(is.na(arriba_plot$color[i]), "grey", arriba_plot$color[i]),
    border = NA
  )
}

# Prepare data for labeling genes at fusion breakpoints
tmp_df1 <- arriba_plot[, c("chr1", "pos1", "pos1", "X.gene1")]
colnames(tmp_df1) <- c("chr", "start", "end", "gene")

tmp_df2 <- arriba_plot[, c("chr2", "pos2", "pos2", "gene2")]
colnames(tmp_df2) <- c("chr", "start", "end", "gene")

# Combine both sides
tmp_df <- rbind(tmp_df1, tmp_df2)

# Add labels if there are any rows
if (nrow(tmp_df) > 0) {
  circos.genomicLabels(
    tmp_df,
    labels.column = "gene",
    side = "inside",
    connection_height = 0.1
  )
}

legend(
  "bottomleft",
  legend = names(type_colors),
  col = type_colors,
  lty = 1,
  lwd = 5,
  cex = 0.8
)
dev.off()


#######################################
# Logging close
#######################################
sink(type = "message")
sink()
close(log_con)
