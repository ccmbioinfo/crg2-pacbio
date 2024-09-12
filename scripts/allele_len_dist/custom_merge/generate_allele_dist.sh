#!/bin/bash

mkdir jobs
mkdir logs

for i in {1..22} X Y
do
	script=jobs/allele_dist_chr${i}.sh
	echo "#!/bin/bash" > $script
	echo "source activate  /hpf/largeprojects/ccmbio/mcouse/conda_envs/notebook_env/" >> $script
	echo "python3 allele_distribution_db.py --alleles_path C4R.alleles.chr${i}.db" >> $script
	sbatch -J allele_dist -o logs/allele_dist_chr${i}.log -e logs/allele_dist_chr${i}.log --nodes=1 --cpus-per-task=1 --mem=40G  --time=100:00:00 $script
done
