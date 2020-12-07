# NOTATES development

## Major Changes
- Changed germline variant filtering (for reporting) to (1) NOT "benign"/"likely benign" in ClinVar (2) MAF < 1% (in gnomAD, ExAC, TGP and ESP) and (3) Non-syn. impact and (4) NOT in FLAGs

## Minor changes and bug fixes

***

# NOTATES v1.5.0

## Major Changes
- added structural variant calling with DELLY (annotation with StructuralVariantAnnotation and VariantAnnotation)
- added "GOAL" exome capture kit information
- updated germline report filtering to report "pathogenic"/"likely pathogenic" variants

## Minor changes and bug fixes
- fixed behavior when no loh is found
- fixed issue in reading enrichment file
- fixed issue in summary report. When there is no gene-level SCNA, message is displayed instead of table
- fixed issue in summary report. When there is no somatic short variants, message is displayed instead of table

***

# NOTATES v1.4.0

## Major Changes
- created environment yaml files (conda) for portability
- changed reference genome to hg38
- changed annotation tool from oncotator to funcotator
- added installation instructions

## Minor changes and bug fixes

***

# NOTATES v1.3.0

## Major Changes
- Changed germline variant filtering (for reporting) to (1) NOT "benign"/"likely benign" in ClinVar (2) MAF < 1% and (3) Non-syn. impact and (4) NOT in FLAGs
- Changed somatic variant filtering (Mutect2) to the new recommendation by GATK
- Added summary report
- Removed usage of compound filtering expressions in germline VCF filtering (For such expressions, if a record is null for or missing a particular annotation in the expression, the tool negates the entire compound expression and so automatically passes the variant record even if it fails on one of the expressions.)

## Minor changes and bug fixes
- germline filter list altered so that up-to-date DDR genes are used
- minor changes in report text

***

# NOTATES v1.2

## Major Changes
- replaced "Fanconi Anemia" genes with "DDR genes" in germline variant reporting
- Updated table of DNA Damage Repair genes
- Updated ClinVar pathogenic variants for germline variants

## Minor changes and bug fixes
- in `Report.Rmd`, changed some subsections' titles
- in `Report.Rmd`, turned all subsection titles to title-case
- fixed issue in enrichment visualization

***

# NOTATES v1.1.4

- updated ClinVar pathogenic
- fixed minor issues
- replaced altered prop. with WGII

***

# NOTATES v1.1.3

minor changes:

- added VAF filter (>0.05) for pathfindR
- added VAF filter (>0.05) for deconstructsigs
- added the WGII as a SCNA metric
- improved logR cut-off for SCNA filter (logR < .25 | >.25)
- other minor changes

***

# NOTATES v1.1.2

fix inflater/deflator issue in MacOSs

***

# NOTATES v1.1.1

fix bug in pathfindR enrichment script

***

# NOTATES v1.1

Fix exome intervals

***

# NOTATES v1.0.1

minor changes/additions and bug fixes

***

# NOTATES v1.0

Reads-to-filtered-variants pipeline for tumor + normal WES analysis of neurooncological cases