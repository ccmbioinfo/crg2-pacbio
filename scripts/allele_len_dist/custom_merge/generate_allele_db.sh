#!/bin/bash
#SBATCH --job-name=allele_db
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --output=%x-%j.out


source activate  /hpf/largeprojects/ccmbio/mcouse/conda_envs/notebook_env/

python3 generate_allele_db.py \
        --vcf_path ../../data/C4R_all_TRGTv0.5.0 \
        --output_file C4R.alleles.db.gz

zcat C4R.alleles.db.gz | sort -T . -k 1,1 | gzip > C4R.alleles.sorted.db.gz
for i in {1..22} X Y;
do
        zcat C4R.alleles.sorted.db.gz | grep chr${i} > C4R.alleles.chr${i}.db
done
