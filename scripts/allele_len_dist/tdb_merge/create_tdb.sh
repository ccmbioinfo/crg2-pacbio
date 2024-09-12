# Generate tdb for each sample (i.e. each TRGT VCF)

# CMH vcfs: /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT
# TCAG vcfs: /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/TCAG_TRGT
vcf_dir=$1

# create list of VCFs
ls ${vcf_dir}/vcf/*.vcf.gz > vcf_paths.txt

mkdir -p ${vcf_dir}/tdbs/
mkdir -p ${vcf_dir}/jobs/

cat vcf_paths.txt | while read in_vcf
do
    sample_name=$(basename $in_vcf)
    sample_name=${sample_name%.vcf.gz}
    
    script_name=${vcf_dir}/jobs/per_sample.${sample_name}.sh
    echo '#!/bin/bash' > $script_name
    # This is just the command from Step 1
    echo "source activate /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/" >> ${script_name} 
    echo "tdb create -o ${vcf_dir}/tdbs/${sample_name}.tdb $in_vcf" >> ${script_name}
done

mkdir -p ${vcf_dir}/logs/
for i in ${vcf_dir}/jobs/per_sample.*
do 
    name=$(basename $i) 
    sbatch -J $name -o ${vcf_dir}/logs/${name}.tdb.create.log -e ${vcf_dir}/logs/${name}.tdb.create.log --nodes=1 --cpus-per-task=1 --mem=8G $i
done


