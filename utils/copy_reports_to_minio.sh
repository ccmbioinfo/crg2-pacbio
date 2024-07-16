#!/bin/bash

family_dir=$1
mc=/hpf/largeprojects/ccmbio/ccmmarvin_shared/data_transfers/mc

if [ -z $family_dir ]; then
        echo 'Must specify results path, exiting'
        exit
fi

family=`echo $family_dir | cut -d '/' -f1`

$mc mb minio/results-c4r/long_read_seq/reports_2024/${family}

# SNVs and indels
$mc cp ${family}/*/small_variants/coding/${family}/*wes* minio/results-c4r/long_read_seq/reports_2024/${family}

# SVs
$mc cp ${family}/*/sv/${family}.pbsv.2*.csv minio/results-c4r/long_read_seq/reports_2024/${family}

# pathogenic REs
$mc cp ${family}/*/pathogenic_repeats/${family}.known.path.str.loci.2*.csv minio/results-c4r/long_read_seq/reports_2024/${family}

# TCAG CNVs
$mc cp ${family}/*/cnv/${family}*cnv.csv  minio/results-c4r/long_read_seq/reports_2024/${family}

# TO DO: genome-wide repeats
