#!/bin/bash

analyses=$1 # path to Stager analyses_report.csv
DATA_DIR="/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/data"
ANALYSIS_DIR="/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/analyses/"
CRG2_PACBIO=~/crg2-pacbio
MAPPING="/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/family_sample_mapping.txt"

module load python/3.7.1

if [ -z $analyses ]; then
        echo 'Must specify path to Stager analyses_report.csv, exiting'
        exit
fi

while read -r line
do
    family=`echo $line | cut -d] -f2 |  tr '"' ' ' | tr "'[]" ' ' | cut -d, -f2 | sed 's/ //g'`
    analysis=`echo $line | grep -oE '\b[0-9]{5}\b' | sed 's/ //g'`
    sample1=`echo $line | cut -d, -f3 |  tr '"' ' ' | tr "'[]" ' ' | sed 's/ //g'`
    echo $family 

    if [ "$family" != "analysis_state" ]; then
        declare -a crams=()
        for sample in `grep $family $MAPPING | awk '{print $2}'`
        do
                echo $sample
                cram=`ls ${DATA_DIR}/*${sample}/*${sample}*haplotagged.bam`
                crams+=($cram)
        done

        echo $crams

        # create analysis directory for family
        if [ ! -d ${ANALYSIS_DIR}/${family} ]; then
                mkdir -p ${ANALYSIS_DIR}/${family}/${analysis}
        elif [ ! -d ${family}/${analysis} ]; then
                mkdir ${ANALYSIS_DIR}/${family}/${analysis}
        fi

        # get HPO terms and pedigree
        # HPO terms will be derived from the proband no matter which family member is specified
        #python3 ~/crg2/scripts//get_HPO_pedigree.py \
        #        -i ${family}_${sample1} \
        #        -c /home/ccmmarvin/crg2/credentials.csv
        
        # copy pipeline files to analysis directory and modify to add family, HPO, and pedigree 
        cp ${CRG2_PACBIO}/config.yaml ${CRG2_PACBIO}/crg2-pacbio.sh ${CRG2_PACBIO}/slurm_profile/slurm-config.yaml ${ANALYSIS_DIR}/${family}/${analysis}
        sed -i "s/NA12878/$family/" ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        sed -i ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        HPO=`echo /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/annotate_SV/HPO/${family}*with_genes.txt`
        pedigree=`echo /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/code/scripts/misc/${family}*`
        echo $HPO $pedigree
        sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+"  ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        sed -i "s+ped: \"\"+ped: \"${pedigree}\"+"  ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        # add targets to job submission script
        #sed -i "s+{SLURM}+{SLURM} -p sv/${family}.pbsv.csv small_variants/coding/${family} pathogenic_repeats/${family}.known.path.str.loci.csv small_variants/panel/${family} small_variants/panel-flank/${family}   pathogenic_repeats/${family}.known.path.str.loci.csv repeat_outliers/${family}.repeat.outliers.annotated.csv TRGT_denovo/${family}.TRGT.denovo.annotated.csv+g" ${ANALYSIS_DIR}/${family}/${analysis}/crg2-pacbio.sh


        # create samples.tsv and copy per-sample TRGT VCFs 
        echo -e "sample\tBAM\tcase_or_control" > ${ANALYSIS_DIR}/${family}/${analysis}/samples.tsv
        echo -e "family\tplatform\tsmall_variant_vcf\ttrgt_vcf_dir\tpbsv_vcf\ttrgt_pathogenic_vcf_dir" > ${ANALYSIS_DIR}/${family}/${analysis}/units.tsv
        for cram in ${crams[@]}
        do
                sample=`basename $cram | cut -d '.' -f1 | cut -d '_' -f2`
                sample=`echo $sample | sed 's/-ready//'`
                echo -e "$sample\t$cram" >> ${ANALYSIS_DIR}/${family}/${analysis}/samples.tsv
		if [ ! -d ${ANALYSIS_DIR}/${family}/${analysis}/trgt ]; then 
			mkdir ${ANALYSIS_DIR}/${family}/${analysis}/trgt
		fi
		trgt=`echo /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/TRGT/genotype_GIAB_repeats_trgtv0.5.0/*${sample}.vcf.gz`
		cp $trgt ${ANALYSIS_DIR}/${family}/${analysis}/trgt
        done

        # create units.tsv
        small_var_vcf=`echo ${DATA_DIR}/cohorts/*${family}/*${family}.GRCh38.deepvariant.glnexus.phased.vcf.gz`
        sv_vcf=`echo ${DATA_DIR}/cohorts/*${family}/*${family}.GRCh38.pbsv.vcf.gz`
        echo -e "$family\tPACBIO\t${small_var_vcf}\ttrgt\t${sv_vcf}\t" >> ${ANALYSIS_DIR}/${family}/${analysis}/units.tsv

        cd $dir
    fi
done<$analyses
