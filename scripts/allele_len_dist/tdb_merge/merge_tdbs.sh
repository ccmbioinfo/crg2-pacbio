# Lets assume there are 100 tdbs and we want to make 10 'level1' merges, each holding 10 tdbs

# CMH tdbs: /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT/tdbs
# TCAG tdbs: /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/TCAG_TRGT/tdbs
# collect tdb paths
tdb_dir=$1
echo ${tdb_dir}/*tdb | sed 's/ /\n/g' > ${tdb_dir}/tdb_paths.txt
split -l 10 ${tdb_dir}/tdb_paths.txt ${tdb_dir}/level1.

# makes e.g. level1.aa, level1.ab, etc

mkdir -p ${tdb_dir}/level1_tdbs/jobs
mkdir -p ${tdb_dir}/level1_tdbs/logs
mkdir -p ${tdb_dir}/level1_tdbs/tdbs

for part_list in ${tdb_dir}/level1.*
do
    partname=`basename $part_list`
    script_name=${tdb_dir}/level1_tdbs/jobs/level1_${partname}.sh
    echo '#!/bin/bash' > $script_name
    echo "source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/" >> $script_name 
    echo "cat $part_list | xargs tdb merge --mem 20 --threads 8 -o ${tdb_dir}/level1_tdbs/tdbs/${partname}.tdb" >> $script_name
    merge_job=`sbatch -J level1 -o ${tdb_dir}/logs/level1_${partname}.log -e ${tdb_dir}/level1_tdbs/logs/level1_${partname}.log --nodes=1 --cpus-per-task=8 --mem=20G $script_name | cut -d' ' -f4`
    echo $merge_job
done

### Step 4: make the final merge

# This does a similar thing to Step 3, just with the `tdb_paths.txt` being filled with 'level1' tdbs
# mkdir -p ${tdb_dir}/level2_tdbs/jobs
# mkdir -p ${tdb_dir}/level2_tdbs/logs
# mkdir -p ${tdb_dir}/level2_tdbs/tdb
# script_name=${tdb_dir}/level2_tdbs/jobs/level2.sh
# echo '#!/bin/bash' > $script_name
# echo "source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/" >> ${script_name} 
# echo "tdb merge --mem 40 --threads 8 -o ${tdb_dir}/level2_tdbs//final_merge.tdb ${tdb_dir}/level1_tdbs/*.tdb" >> $script_name
# sbatch -J level2 -o ${tdb_dir}/level2_tdbs/logs/level2.log -e ${tdb_dir}/level2_tdbs/logs/level2.log --nodes=1 --cpus-per-task=8 --mem=20G --dependency=afterok:$merge_job  $script_name
