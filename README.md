# <img src="Scripts/NOT_logo.png" align="left" height=120/> NOTATES : NOT - Alterations in Tumor Exome Sequencing
NeuroOncology Technologies(NOT) aims to analyze Whole Exome Sequencing data of blood-matched tumor samples that belong to brain tumor patients. For the moment, the analysis is focused on gliomas and medulloblastomas (MBs). Therefore, glioma/MB-specific genomic alterations that are expected to affect diagnosis, prognosis and treatment response are assessed in each case, utilizing manually-curated databases of genes and genomic alterations that are related to brain tumor biology.

# Overview of the Pipeline
The pipeline consists of several bash and R scripts, and a wrapper bash script.

The combined germline variation, somatic Single Nucleotide Variant(SNV) and Insertion/Deletion(InDel) discovery pipeline was created based on the GATK Best Practices workflow. Germline SNV and indels are called using GATK HaplotypeCaller. Somatic SNV and indels are called using GATK MuTect2. Both germline and somatic variants are annotated via Funcotator.

Somatic copy number alterations (SCNA) and loss of heterozygosity (LOH) events are called using ExomeCNV (Sathirapongsasuti JF, Lee H, Horst BA, et al. Exome sequencing-based copy-number variation and loss of heterozygosity detection: ExomeCNV. Bioinformatics. 2011;27(19):2648-54).

Clonal/subclonal copy number aberrations are estimated using THetA2 (Oesper L, Satas G, Raphael BJ. Quantifying tumor heterogeneity in whole-genome and whole-exome sequencing data. Bioinformatics. 2014;30(24):3532-40).

MSI status is predicted using MSIpred (Wang C, Liang C. MSIpred: a python package for tumor microsatellite instability classification from tumor mutation annotation data using a support vector machine. Sci Rep. 2018;8(1):17546).

Mutational signatures are estimated using DeconstructSigs (Rosenthal R, Mcgranahan N, Herrero J, Taylor BS, Swanton C. DeconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. Genome Biol. 2016;17:31).

Pathway enrichment analyses are performed using [pathfindR](https://github.com/egeulgen/pathfindR) (Ulgen E, Ozisik O, Sezerman OU. 2019. pathfindR: An R Package for Comprehensive Identification of Enriched Pathways in Omics Data Through Active Subnetworks. Front. Genet. https://doi.org/10.3389/fgene.2019.00858)

Variant reporting and prioritization of detected alterations are achieved via custom scripts and curated resources.

# Installation Instructions

See [here](installation_instructions.md) for installation instructions.
