#!/bin/bash
#SBATCH --job-name=crg2-pacbio
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out

SF="/hpf/largeprojects/ccmbio/mcouse/pacbio_report_dev/tools/crg2-pacbio/Snakefile"
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake"
SLURM="/hpf/largeprojects/ccmbio/mcouse/pacbio_report_dev/tools/crg2-pacbio/slurm_profile/"
CONFIG="config.yaml"

module purge

source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --configfile ${CONFIG} --profile ${SLURM}