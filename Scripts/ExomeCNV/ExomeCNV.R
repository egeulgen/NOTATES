############################################################
#                NeuroOncology Technologies                #
#             Whole-Exome Sequencing Pipeline              #
#                    ExomeCNV Analysis                     #
#                   Ege Ulgen, Oct 2016                    #
############################################################

require(ExomeCNV)
require(DNAcopy)

# set workdir to currentdir/ExomeCNV
setwd("./ExomeCNV/")

# LOH Calling -------------------------------------------------------------
load("./baf/BAF_data.Rdata")

### Calling LOH on each heterozygous position
eLOH <- LOH.analyze(normal = normal_BAF, tumor = tumor_BAF, alpha=0.05, method="two.sample.fisher")

###Combine multiple positions into LOH segments
the.loh <- multi.LOH.analyze(normal = normal_BAF, tumor = tumor_BAF, 
                             all.loh.ls = list(eLOH), 
                             test.alpha=0.001, method="variance.f", 
                             sdundo=c(0,0), alpha=c(0.05,0.01))

pdf("LOH.pdf", width = 15, height = 8)
do.plot.loh(the.loh, normal_BAF, tumor_BAF, "two.sample.fisher", plot.style="baf")
dev.off()

LOH.regions <- the.loh[the.loh$LOH, ]
write.csv(LOH.regions, file = "LOH_regions.csv", row.names = F)

tumor_loh_sites <- expand.loh(the.loh, tumor_BAF)
tumor_loh_sites <- tumor_loh_sites[tumor_loh_sites$LOH, ]

normal_loh_sites <- expand.loh(the.loh, normal_BAF)
normal_loh_sites <- normal_loh_sites[normal_loh_sites$LOH, ]

LOH.sites <- merge(normal_loh_sites, tumor_loh_sites, by=c("chr","position"), sort = F)
cnames <- colnames(LOH.sites)
cnames <- sub(".x", "_normal", cnames)
cnames <- sub(".y", "_tumor", cnames)
colnames(LOH.sites) <- cnames
LOH.sites <- LOH.sites[,!(grepl("LOH", colnames(LOH.sites)))]
write.csv(LOH.sites, "LOH_sites.csv", row.names = F)

# CNV Calling -------------------------------------------------------------
chr.list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

###### Load in coverage files
read.coverage.gatk.fix <- function(file){
  gatk = read.table(file, header = TRUE)
  gatk <- gatk[grepl("-", gatk$Target), ]
  
  chrpos = matrix(unlist(strsplit(as.character(gatk$Target), 
                                  ":")), ncol = 2, byrow = TRUE)
  chr = factor(chrpos[, 1])
  pos = matrix(as.integer(unlist(strsplit(chrpos[, 2], "-"))), 
               ncol = 2, byrow = TRUE)
  start = pos[, 1]
  end = pos[, 2]
  
  return(data.frame(probe = gatk$Target, chr = chr, probe_start = start, 
                    probe_end = end, targeted.base = end - start + 1, sequenced.base = NA, 
                    coverage = as.numeric(gatk$total_coverage), average.coverage = as.numeric(gatk$average_coverage), 
                    base.with..10.coverage = NA))
}
normal <- read.coverage.gatk.fix("./DepthOfCoverage/normal.coverage.sample_interval_summary")
tumor <- read.coverage.gatk.fix("./DepthOfCoverage/tumor.coverage.sample_interval_summary")

###### Calculate log coverage ratio
patient.logR <- calculate.logR(normal, tumor)

# admix_rate <- 1 - 2*mean(sapply(LOH.regions$tumor.baf/LOH.regions$tumor.coverage, function(x) abs(x - 0.5)))
# print(admix_rate)
admix_rate <- 1 - 2*mean(sapply(LOH.sites$baf_tumor/LOH.sites$coverage_tumor, function(x) abs(x - 0.5)))
print(admix_rate)


###### Call CNV for each exon
### Call CNV on each exon (using classify.eCNV), one chromosome at a time. We recommend high min.spec (0.9999) 
### and option="spec" to be conservative against false positive. This is because whatever is called 
### at exon level will persist through merging step
patient.eCNV <- c()
for (i in 1:length(chr.list)) {
  idx <- (normal$chr == chr.list[i])
  ecnv <- classify.eCNV(normal = normal[idx, ], tumor = tumor[idx, ], logR=patient.logR[idx], 
                        min.spec = 0.9999, min.sens = 0.9999, option="spec", admix = admix_rate, read.len = 76)
  patient.eCNV <- rbind(patient.eCNV, ecnv)
}

###### Combine exonic CNV into larger segments
###Here, we use lower min.spec and min.sens and option="auc" to be less conservative and allow for more discovery.
patient.cnv <- multi.CNV.analyze(normal, tumor, logR = patient.logR, all.cnv.ls = list(patient.eCNV), 
                                 min.spec = 0.99, min.sens = 0.99, option = "auc", 
                                 coverage.cutoff = 5, admix = admix_rate, read.len = 76)

###### plot the results and export outputs.
pdf("CNV.pdf", width = 15, height = 8)
do.plot.eCNV(patient.cnv, style = "bp", bg.cnv = patient.eCNV, line.plot = T)
dev.off()

write.output(patient.eCNV, patient.cnv, "CNV")
