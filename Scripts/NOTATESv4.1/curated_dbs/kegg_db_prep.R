pw_ids <- c("hsa04110", "hsa04150", "hsa05200", 
            "hsa04310", "hsa04340")
names(pw_ids) <- c("Cell Cycle", "mTOR Signalling", "Pathways in Cancer",
                   "WNT Signalling", "Hedgehog Signalling")

final <- c()
for (pw_id in pw_ids) {
  pw <- KEGGREST::keggGet(pw_id)
  pw <- pw[[1]]$GENE[c(FALSE, TRUE)] ## get gene symbols, not descriptions
  pw <- sub(";.+", "", pw) ## discard any remaining description
  pw <- pw[grepl("^[a-zA-Z0-9_-]*$", pw)] ## remove mistaken lines that cannot be gene symbols
  
  final <- rbind(final, data.frame(Gene = pw, pathway = names(pw_ids)[pw_ids == pw_id]))
}

write.csv(final, "~/Downloads/kegg_pws.csv", row.names = FALSE)
