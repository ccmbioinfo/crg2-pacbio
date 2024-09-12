#!/bin/bash
#SBATCH --job-name=parse_TRGT
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80G
#SBATCH --output=%x-%j.out

input_VCF=$1
output_TSV=$2

module load python/3.7.1
echo -e "Starting run at: `date`"
python3 parse_TRGT_multi-sample_VCF.py --input_vcf $input_VCF --output_tsv $output_TSV
echo -e "Job finished with exit code $? at: `date`"