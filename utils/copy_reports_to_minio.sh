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

# panel report
panel=${family}/*/small_variants/panel/*/*wgs*
panel_date=`basename $panel | cut -d'.' -f3`
$mc cp $panel minio/results-c4r/long_read_seq/reports_2024/${family}/${family}.wgs.panel.${panel_date}.csv

# panel flank report
panelflank=${family}/*/small_variants/panel-flank/*/*wgs*
panel_date=`basename $panelflank | cut -d'.' -f3`
$mc cp $panelflank minio/results-c4r/long_read_seq/reports_2024/${family}/${family}.wgs.panel-flank.${panel_date}.csv

# SVs
$mc cp ${family}/*/sv/${family}.pbsv.2*.csv minio/results-c4r/long_read_seq/reports_2024/${family}

# pathogenic REs
$mc cp ${family}/*/pathogenic_repeats/${family}.known.path.str.loci.2*.csv minio/results-c4r/long_read_seq/reports_2024/${family}

# repeat expansion outliers
$mc cp ${family}/*/repeat_outliers/${family}.repeat.outliers.annotated.2*.csv minio/results-c4r/long_read_seq/reports_2024/${family}

# TRGT de novo tandem repeats
$mc cp ${family}/*/TRGT_denovo/${family}*TRGT.denovo.annotated.2*.csv minio/results-c4r/long_read_seq/reports_2024/${family}

# TCAG CNVs
$mc cp ${family}/*/cnv/${family}*cnv*csv  minio/results-c4r/long_read_seq/reports_2024/${family}

