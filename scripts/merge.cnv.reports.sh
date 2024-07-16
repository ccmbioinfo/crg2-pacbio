#!/bin/bash
#SBATCH --job-name=cnv
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out

module load python/3.7.1

PROJECT=$1
TODAY=`date +%Y-%m-%d`
OUT="${PROJECT}.${TODAY}.cnv.tsv"

if [[ -z "$PROJECT" ]]; then
	PROJECT="PROJECT";
fi

python3 ~/crg/merge.cnv.reports.py -i $(ls *.tsv | tr '\n' ' ') -o ${OUT}
python3 ~/crg/cap_report.py ${OUT}
