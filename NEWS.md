# NOTATES v1.3.0.9000
# to be released as v1.3.0
## Major Changes
- Changed germline filtering to (1) NOT "benign"/"likely benign" in ClinVar (2) MAF < 1% and (3) Non-syn. impact and (4) NOT in FLAGs

## Minor changes and bug fixes
- germline filter list altered so that up-to-date DDR genes are used

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