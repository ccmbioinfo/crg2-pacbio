#!/bin/bash

analyses=$1 # path to Stager analyses_report.csv
TCAG_SHARED="/hpf/largeprojects/ccmbio/ccmmarvin_tcag_shared/PacBio"
ANALYSIS_DIR="/hpf/largeprojects/ccmbio/mcouse/pacbio_report_dev/results/test_setup"
CRG2_PACBIO=~/crg2-pacbio

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
        family_=`echo $family _ | sed 's/ //'`
        crams=`find $TCAG_SHARED -name '*haplotagged.bam' | grep -E $family_`
        echo $crams

        # create analysis directory for family
        if [ ! -d ${ANALYSIS_DIR}/${family} ]; then
                mkdir -p ${ANALYSIS_DIR}/${family}/${analysis}
        elif [ ! -d ${family}/${analysis} ]; then
                mkdir ${ANALYSIS_DIR}/${family}/${analysis}
        fi

        # get HPO terms and pedigree
        # HPO terms will be derived from the proband no matter which family member is specified
        python3 ~/crg2/scripts//get_HPO_pedigree.py \
                -i ${family}_${sample1} \
                -c /home/ccmmarvin/crg2/credentials.csv
        
        # copy pipeline files to analysis directory and modify to add family, HPO, and pedigree 
        cp ${CRG2_PACBIO}/config.yaml ${CRG2_PACBIO}/crg2-pacbio.sh ${CRG2_PACBIO}/slurm_profile/slurm-config.yaml ${ANALYSIS_DIR}/${family}/${analysis}
        sed -i "s/NA12878/$family/" ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        sed -i ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        HPO=`echo ~/gene_data/HPO/${family}*`
        pedigree=`echo ~/gene_data/pedigrees/${family}*`
        echo $HPO $pedigree
        sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+"  ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        sed -i "s+ped: \"\"+ped: \"${pedigree}\"+"  ${ANALYSIS_DIR}/${family}/${analysis}/config.yaml
        # add targets to job submission script
        sed -i "s+{SLURM}+{SLURM} -p sv/${family}.pbsv.csv small_variants/coding/${family} pathogenic_repeats/${family}.known.path.str.loci.csv+g" ${ANALYSIS_DIR}/${family}/${analysis}/crg2-pacbio.sh


        # create samples.tsv 
        echo -e "sample\tBAM" > ${ANALYSIS_DIR}/${family}/${analysis}/samples.tsv
        echo -e "family\tplatform\tsmall_variant_vcf\ttrgt_vcf_dir\tpbsv_vcf\ttrgt_pathogenic_vcf_dir" > ${ANALYSIS_DIR}/${family}/${analysis}/units.tsv
        for cram in ${crams[@]}
        do
                echo $cram
                sample=`basename $cram | cut -d '.' -f1 | cut -d '_' -f2`
                sample=`echo $sample | sed 's/-ready//'`
                echo -e "$sample\t$cram" >> ${ANALYSIS_DIR}/${family}/${analysis}/samples.tsv
        done

        # create units.tsv
        small_var_vcf=`echo ${TCAG_SHARED}/*/${family}*joint/${family}.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz`
        if [ ! -f $small_var_vcf ]; then # singleton
                small_var_vcf=`echo ${TCAG_SHARED}/*/${family}*/${family}*GRCh38.deepvariant.phased.vcf.gz`
        fi
        sv_vcf=`echo ${TCAG_SHARED}/*/${family}*joint/${family}.joint.GRCh38.pbsv.phased.vcf.gz`
        if [ ! -f $sv_vcf ]; then # singleton
                sv_vcf=`echo ${TCAG_SHARED}/*/${family}*/${family}*GRCh38.pbsv.phased.vcf.gz`
        fi 
        echo -e "$family\tPACBIO\t${small_var_vcf}\t\t${sv_vcf}\t" >> ${ANALYSIS_DIR}/${family}/${analysis}/units.tsv

        # copy TCAG annotated CNV files for merging
        mkdir ${ANALYSIS_DIR}/${family}/${analysis}/cnv
        cp /hpf/largeprojects/ccmbio/ccmmarvin_tcag_shared/PacBio/*/${family}*/annotation_outputs/*cnv.tagged.tsv ${ANALYSIS_DIR}/${family}/${analysis}/cnv
        dir=`pwd`
        cd ${ANALYSIS_DIR}/${family}/${analysis}/cnv
        sbatch ${CRG2_PACBIO}/scripts/merge.cnv.reports.sh $family
        cd $dir
    fi
done<$analyses
