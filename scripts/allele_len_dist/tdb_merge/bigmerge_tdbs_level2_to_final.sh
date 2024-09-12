#!/bin/bash

tdb_dir=$1
### Make the final merge (if only two merges required, e.g. for TCAG. CMH requires three merges)

mkdir -p ${tdb_dir}/level2_tdbs/jobs
mkdir -p ${tdb_dir}/level2_tdbs/logs
mkdir -p ${tdb_dir}/level2_tdbs/tdb
script_name=${tdb_dir}/level2_tdbs/jobs/level2.sh
echo '#!/bin/bash' > $script_name
echo "source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/" >> ${script_name} 
echo "tdb merge --mem 60 --threads 8 -o ${tdb_dir}/level2_tdbs//final_merge.tdb ${tdb_dir}/*.tdb" >> $script_name
sbatch -J level2 -o ${tdb_dir}/level2_tdbs/logs/level2.log -e ${tdb_dir}/level2_tdbs/logs/level2.log --nodes=1 --cpus-per-task=8 --mem=60G  $script_name
