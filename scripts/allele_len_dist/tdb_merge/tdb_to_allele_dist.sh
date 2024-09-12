#!/bin/bash
#SBATCH --job-name=tdb-to-db
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300G
#SBATCH --output=%x-%j.out

module load anaconda
source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb
python3 tdb_to_allele_dist.py 



