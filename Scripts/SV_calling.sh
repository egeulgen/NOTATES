mkdir DELLY

### Create samples tsv
TUMOR_ID=$(samtools view -H tumor.final.bam | grep "^@RG" | cut -f5 -d$'\t' | sed -e "s/^SM://")
NORMAL_ID=$(samtools view -H normal.final.bam | grep "^@RG" | cut -f5 -d$'\t' | sed -e "s/^SM://")
echo $TUMOR_ID$'\t'tumor > DELLY/samples.tsv
echo $NORMAL_ID$'\t'control >> DELLY/samples.tsv

### Initial SV calling
delly call -x $DELLY_exclude -o DELLY/tumor_sv.bcf -g $genome tumor.final.bam normal.final.bam
### Pre-filter
delly filter -f somatic -o DELLY/tumor_sv_pre.bcf -s DELLY/samples.tsv DELLY/tumor_sv.bcf
### Second SV calling
delly call -g $genome -v DELLY/tumor_sv_pre.bcf -o DELLY/geno.bcf -x $DELLY_exclude tumor.final.bam normal.final.bam
### Post-filter
delly filter -f somatic -o DELLY/somatic_SV.bcf -s DELLY/samples.tsv DELLY/geno.bcf

### BCF to VCF
bcftools view DELLY/somatic_SV.bcf > DELLY/somatic_SV.vcf
