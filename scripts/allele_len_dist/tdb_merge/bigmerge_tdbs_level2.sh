#!/bin/bash
# merge single-sample tdbs into 10-sample tdbs

# tdb_dir=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT/tdbs
tdb_dir=$1

export TMPDIR=.


echo ${tdb_dir}/level1_tdbs/tdbs/*tdb | sed 's/ /\n/g' > ${tdb_dir}/level1_tdbs/tdb_paths.txt
# makes e.g. level1.aa, level1.ab, etc
split -l 10 ${tdb_dir}/level1_tdbs/tdb_paths.txt ${tdb_dir}/level1_tdbs/level2.

mkdir -p ${tdb_dir}/level2_tdbs/jobs
mkdir -p ${tdb_dir}/level2_tdbs/logs
mkdir -p ${tdb_dir}/level2_tdbs/tdbs

for part_list in ${tdb_dir}/level1_tdbs/level2.*
do
    partname=`basename $part_list`
    script_name=${tdb_dir}/level2_tdbs/jobs/level2_${partname}.sh
    echo '#!/bin/bash' > $script_name
    echo 'export TMPDIR=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/code/scripts/' >> $script_name
    echo "source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/" >> $script_name
    echo "cat $part_list | xargs tdb bigmerge --mem 80 --threads 12 -o ${tdb_dir}/level2_tdbs/tdbs/${partname}.tdb" >> $script_name
    sbatch -J level2 -o ${tdb_dir}/level2_tdbs/logs/level2_${partname}.log -e ${tdb_dir}/level2_tdbs/logs/level2_${partname}.log --nodes=1 --cpus-per-task=12 --mem=80G $script_name
done
