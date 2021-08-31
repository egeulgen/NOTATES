# 1. Install conda

See instructions on https://docs.conda.io/en/latest/miniconda.html

# 2. Create conda environments

```bash
conda env create -f NOTATES_main_env.yml
conda env create -f NOTATES_R_env.yml
conda env create -f NOTATES_python_env.yml
```

# 3. Activate NOTATES main environment and register GATK3

GATK3 jar file can be found here: https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk

``` bash
conda activate NOTATES_main
gatk3-register /path/to/GenomeAnalysisTK.jar 
conda deactivate
```

# 4. Activate NOTATES R environment and install non-conda R packages

``` bash
conda activate NOTATES_R
R
```

``` R
install.packages("https://cran.r-project.org/src/contrib/Archive/ExomeCNV/ExomeCNV_1.4.tar.gz", 
                  repos = NULL, type = "source", method = "libcurl")
devtools::install_github("raerose01/deconstructSigs")
devtools::install_github("egeulgen/pathfindR") # will take a long time
```

# 5. Install THetA2

Activate the NOTATES_python environment, obtain the most recent version and install:

``` bash
conda activate NOTATES_python
git clone https://github.com/raphael-group/THetA.git
cd THetA
./install
conda deactivate
```

# 6. Install MSIpred

Activate the NOTATES_python environment, obtain the most recent version and install:

``` bash
conda activate NOTATES_python
git clone https://github.com/wangc29/MSIpred.git
cd MSIpred
python setup.py install
conda deactivate
```

# 7. Download necessary resources

The pipeline utilizes the following resources (`path/to/resource` can be altered in the configuration file `notates.config`):

- Reference sequence, prepared for use with BWA and GATK (can be found in the [GATK bundle])
- Exome capture kit probe and target intervals (check with your capture kit provider or your sequencing facility)
- dbSNP VCF (can be found in the [GATK bundle])
- Mills and 1000 Genomes gold standard InDel Sites VCF (can be found in the [GATK bundle])
- 1000 Genomes phase 1 high confidence SNP Sites VCF (can be found in the [GATK bundle])
- `af-only-gnomad` VCF (can be found in [GATK somatic bundle])
- `small_exac_common` VCF (can be found in [GATK somatic bundle])
- Funcotator data sources (see the [tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial))

[GATK bundle]: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
[GATK somatic bundle]: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38

Additionally, the following files should be placed in the designated places:

- Table of Cancer Gene Census genes (to be placed under `Scripts/CGC_latest.csv`) (can be obtained for non-commercial purposes on [https://cancer.sanger.ac.uk/census](https://cancer.sanger.ac.uk/census))
- The table of COSMIC SBS Signatures descriptions (to be placed under `Scripts/cosmic_sbs_details.csv`, required columns are: Signature, Proposed aetiology, Associated mutation classes and signatures, Comments) (can be obtained for non-commercial purposes on [https://cancer.sanger.ac.uk/cosmic/signatures/index.tt](https://cancer.sanger.ac.uk/cosmic/signatures/index.tt))
- Tables of selected KEGG pathways (to be placed under  `Scripts/NOTATES/curated_dbs/important_KEGG_pws_glioma.csv` for glioma and `Scripts/NOTATES/curated_dbs/important_KEGG_pws_MB.csv` for medulloblastoma, with 2 columns: "Gene" and "pathway") (The R script in `Scripts/NOTATES/curated_dbs/kegg_db_prep.R` can be used to obtain the genes under selected KEGG pathways, only for non-commercial purposes)

# 8. Update the configuration file

Update the configuration file `notates.config`

# 9. Run run_me.sh

``` bash
bash run_me.sh Analyis_ID Normal_ID Tumor_ID Capture_kit_name type_of_tumor["glioma" or "MB"] primary?["TRUE" or "FALSE"] tumor_sample_type["LiN2" or "FFPE"]
```
