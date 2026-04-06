1. download reference genome fasta and GTF files. e.g.:

```bash
## H_sapiens
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
```

2. Install the main `conda` environment for the pipeline:

Install the 'notates-rnaseq-pipeline' `conda` environment using the yaml file (from the base directory):

```bash
conda env create --file environment.yaml
```

4. Run the pipeline:

Before beginning make the pipeline wrapper script executable via:

```bash
chmod +x run_notates_rnaseq.py
```

You may now directly run the pipeline wrapper script:

```bash
# activate environment
conda activate notates-rnaseq-pipeline

# run wrapper script (example command below)
./run_notates_rnaseq.py \
	--threads 10 \
	--data-dir data/P1136_B \
	--sample-id NOT-0123-Tumor \
	--reference-genome reference_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--gtf-file reference_data/Homo_sapiens.GRCh38.115.gtf \
	--fusion-blacklist reference_data/blacklist_hg38_GRCh38_v2.5.1.tsv.gz \
	--known-fusions reference_data/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
	--star-index reference_data/STAR_idx \
	--rsem-index reference_data/rsem_idx \
	--output-dir results \
	--create-DAG
```

The FASTQ file names should be as follows:

- `samplename_1.fastq.gz` and `samplename_2.fastq.gz`
		
