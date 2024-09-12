#!/bin/bash
#SBATCH --job-name=trgt_merge
#SBATCH --time=20:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --output=%x-%j.out

# example usage: sbatch TRGT_merge_VCFs.sh /path/to/TRGTv1.1.2/trgt  /path/to/VCFs/ /path/to/reference/human_GRCh38_no_alt_analysis_set.fasta multi_sample_TRGT
export TMPDIR=.
TRGT=$1 # path to TRGT executable, must be version 1.1.0 or higher 
VCF_PATH=$2 # path to TRGT VCFs, e.g. /path/to/VCFs/
echo $VCF_PATH
REF=$3 # path to reference genome fasta
OUTPUT_PREFIX=$4
$TRGT -vv merge --vcf ${VCF_PATH}/*sorted*vcf.gz --genome $REF --output-type z > ${OUTPUT_PREFIX}.vcf.gz

module load bcftools
bcftools sort -O z -o ${OUTPUT_PREFIX}.sorted.vcf.gz ${OUTPUT_PREFIX}.vcf.gz
tabix ${OUTPUT_PREFIX}.sorted.vcf.gz
rm ${OUTPUT_PREFIX}.vcf.gz
