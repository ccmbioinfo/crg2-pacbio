#!/bin/bash

# tdb_dir=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT/tdbs
tdb_dir=$1
# make final merge of level 2 tdbs to tdb containing all samples
export TMPDIR=.

if [ -d /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT/tdbs/level2_tdbs/tdb/final_merge.tdb ]; then
    rm -r /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT/tdbs/level2_tdbs/tdbs/final_merge.tdb
fi

mkdir ${tdb_dir}/final_tdb/

script_name=${tdb_dir}/final_tdb/level3.sh
echo '#!/bin/bash' > $script_name
echo 'export TMPDIR=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/code/scripts/' >> $script_name
echo "source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/" >> $script_name
echo "tdb bigmerge --mem 240 --threads 8 -o ${tdb_dir}/final_tdb/final.tdb ${tdb_dir}/level2_tdbs/tdbs/*.tdb" >> $script_name
sbatch -J level3 -o ${tdb_dir}/final_tdb/level3.log -e ${tdb_dir}/final_tdb/level3.log --nodes=1 --cpus-per-task=8 --mem=240G $script_name

