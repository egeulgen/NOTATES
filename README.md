# NOTATES : NOT - Alterations in Tumor Exome Sequencing
NeuroOncology Technologies(NOT) aims to analyze Whole Exome Sequencing reads of blood-matched tumor samples that belong to 
brain tumor patients. For the moment, the analyses is focused on gliomas. Therefore, glioma-specific genomic alterations that are
expected to affect diagnosis, prognosis and treatment response are assessed in each case, utilizing manually-curated databases of 
genes and genomic alterations concerning glioma biology.
## Overview of the Pipeline
The current somatic SNV and InDel discovery pipeline was created based upon the GATK Best Practices workflow.
### Features to be added
- [x] Update COSMIC, CGC, GATK Bundle, oncotator db, and dbSNP databases, list of "pathogenic" and/or "drugresponse" dbSNP entries
- [x] Update tools GATK, oncotator, samtools, PICARD
- [ ] Update CNV-calling (keep ExomeCNV atm) - for.loh.pl, admix. rate, logR, seg. annot. 
      also add TheTa (n=2) results to final integration
- [ ] Update Germnline parsing
- [ ] Automate QC table creation
- [ ] Create and integrate PoN
- [ ] Integrate ConTest or a similar method for contamination estimation
- [ ] Add a tool for SV discovery
- [ ] HTML report generation
- [ ] Drug targets (determine using hypergeometric testing)
