# XGEN --------------------------------------------------------------------
exome <- read.delim("~/NOTATES/Resources/Data_sources/Exome_Intervals/xGen/xgen-exome-research-panel-probes.bed", header = F)
sum(exome$V3 - exome$V2 + 1) / 10^6
# 51.740863


# Nimblegen ---------------------------------------------------------------
exome <- read.delim("~/NOTATES/Resources/Data_sources/Exome_Intervals/SeqCap_EZ_v2/SeqCap_EZ_Exome_v2_capture_coord.bed", header = F)
sum(exome$V3 - exome$V2 + 1) / 10^6
# 44.23414

# Medexome ---------------------------------------------------------------
exome <- read.delim("~/NOTATES/Resources/Data_sources/Exome_Intervals/MedExome/MedExome_hg19_capture_targets.bed", header = F)
sum(exome$V3 - exome$V2 + 1) / 10^6
# 46.80795