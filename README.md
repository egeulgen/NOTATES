# <img src="Scripts/NOT_logo.png" align="left" height=120/> NOTATES : NOT - Alterations in Tumor Exome Sequencing
NeuroOncology Technologies(NOT) aims to analyze Whole Exome Sequencing reads of blood-matched tumor samples that belong to brain tumor patients. For the moment, the analysis is focused on gliomas and medulloblastomas (MBs). Therefore, glioma/MB-specific genomic alterations that are expected to affect diagnosis, prognosis and treatment response are assessed in each case, utilizing manually-curated databases of genes and genomic alterations that are related to glioma/MB biology.

## Overview of the Pipeline
The pipeline consists of several bash and R scripts, and a wrapper bash script.

The combined germline variation, somatic Single Nucleotide Variant(SNV) and Insertion/Deletion(InDel) discovery pipeline was created based on the GATK Best Practices workflow. Germline SNV and InDels are called using GATK HaplotypeCaller. Somatic SNV and InDels are called using GATK MuTect2. Both germline and somatic variants are annotated via Oncotator.

Somatic copy number alterations (SCNA) are called using ExomeCNV.

Clonal/subclonal copy number aberrations are estimated using THetA2.

Variant reporting and integration of detected alterations are achieved with custom scripts and curated resources.

### Dependencies
The pipeline depends on the following (if not in bin, path/to/tool can be altered in the configuration file configurations.cfg):

Should be located in bin directory:

- [JAVA](https://www.java.com/en/download/manual.jsp)
- [R](https://www.r-project.org)
- [python 2.7.x](https://www.python.org/downloads/)
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [SAMTools](http://samtools.sourceforge.net/)
- [BWA](http://bio-bwa.sourceforge.net/)

For the report, the following should be installed:

- [pandoc](http://pandoc.org/)
- TeX (MacTeX on Mac OS, TeX Live on Linux)

The the `path/to/tool` for the following should be specified in the configuration file:

- [Picard](http://broadinstitute.github.io/picard/index.html)
- [GATK](https://software.broadinstitute.org/gatk/)
- [THetA2](http://compbio.cs.brown.edu/projects/theta/)
- [Oncotator](https://github.com/broadinstitute/oncotator/releases)

### Resources
Note that **only hg19 reference and resources with hg19 coordinates are supported** as this is the only reference genome build compatible with oncotator data sources.

The pipeline utilizes the following resources (`path/to/resource` can be altered in the configuration file _configurations.cfg_):

- Reference sequence, prepared for use with BWA and GATK
- Exome capture kit target intervals (check with your capture kit provider or your sequencing facility)
- dbSNP VCF (can be found in the [GATK bundle])
- Mills and 1000 Genomes gold standard InDel Sites VCF (can be found in the [GATK bundle])
- 1000 Genomes phase 1 InDel Sites VCF (can be found in the [GATK bundle])
- COSMIC VCF - coding and non-coding mutations (follow the instructions on [COSMIC](http://cancer.sanger.ac.uk/cosmic/download), then merge coding and non-coding VCF files and liftover to hg19)
- Oncotator data sources (latest can be dowloaded from [here](https://personal.broadinstitute.org/lichtens/oncobeta/oncotator_v1_ds_Jan262015.tar.gz))

[GATK bundle]: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle

Additionaly, the following files should be placed in the designated places:

- Table of Cancer Gene Census genes (to be placed under `Scripts/CGC_latest.csv`) (can be obtained for non-commercial purposes on [https://cancer.sanger.ac.uk/census](https://cancer.sanger.ac.uk/census))
- The table of COSMIC SBS Signatures descriptions (to be placed under `Scripts/cosmic_sbs_details.csv`, required columns are: Signature, Proposed aetiology, Associated mutation classes and signatures, Comments) (can be obtained for non-commercial purposes on [https://cancer.sanger.ac.uk/cosmic/signatures/index.tt](https://cancer.sanger.ac.uk/cosmic/signatures/index.tt))
- Tables of selected KEGG pathways (to be placed under `Scripts/NOTATES/curated_dbs/important_KEGG_pws.csv` for glioma and `Scripts/NOTATES/curated_dbs/important_KEGG_pws_MB.csv` for medulloblastoma, with 2 columns: "Gene" and "pathway") (The R script in `Scripts/NOTATES/curated_dbs/kegg_db_prep.R` can be used to obtain the genes under selected KEGG pathways, only for non-commercial purposes)
