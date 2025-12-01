rule merge_cnv_reports:
    input:
        get_cnv_dir
    output:
        cnv_report="cnv/{family}.cnv.csv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    log:
        "logs/cnv/merge_cnv_reports_{family}.log"
    conda:
        "../envs/str_sv.yaml"
    shell:
        '''
        TODAY=`date +%Y-%m-%d`
        OUT="{wildcards.family}.cnv.${{TODAY}}.csv"
        python3 {params.crg2_pacbio}/scripts/merge.cnv.reports.py -i $(ls {input}/*.tsv | tr '\n' ' ') -o $OUT
        ln -s $OUT {output.cnv_report}
        '''