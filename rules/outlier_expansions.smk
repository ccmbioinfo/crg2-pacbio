rule merge_trgt_vcf:
    input: config["run"]["units"]
    params:
        trgt = config["tools"]["trgt"],
        genome = config["ref"]["genome"]
    output: "repeat_outliers/{family}.trgt.vcf.gz"
    log: "logs/repeat_outliers/{family}.merge.trgt.vcf.log"
    shell:
        """
        export TMPDIR=`pwd`
        dir=`awk -F "\t" '{{print $4}}' {input} | tail -n 1`
        vcf=`echo $dir/*vcf.gz`
        {params.trgt} -vv merge \
            --vcf $vcf \
            --force-single \
            --genome {params.genome} \
            --output-type z > {output}
        """

rule sort_merged_trgt_vcf:
    input: "repeat_outliers/{family}.trgt.vcf.gz"
    output: "repeat_outliers/{family}.trgt.sorted.vcf.gz"
    log: "logs/repeat_outliers/{family}.trgt.sorted.vcf.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        bcftools sort -O z -o {output} {input}
        tabix {output}
        """

rule find_repeat_outliers:
    input: 
        case_vcf = "repeat_outliers/{family}.trgt.sorted.vcf.gz",
        control_vcf = config["trgt"]["control_alleles"]
    output: "repeat_outliers/{family}.repeat.outliers.tsv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    log:  "logs/repeat_outliers/{family}.repeat.outliers.log"
    resources:
        mem_mb = 100000
    conda: 
        "../envs/str_sv.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/find_repeat_outliers.py --case_vcf {input.case_vcf} \
            --control_vcf {input.control_vcf} \
            --output_file  {output}) > {log} 2>&1 
        """


rule annotate_repeat_outliers:
    input: "repeat_outliers/{family}.repeat.outliers.tsv"
    output: "repeat_outliers/{family}.repeat.outliers.annotated.csv"
    log:  "logs/repeat_outliers/{family}.annotate.repeat.outliers.log"
    params: 
      crg2_pacbio = config["tools"]["crg2_pacbio"],
      genes = config["trgt"]["ensembl"],
      OMIM = config["annotation"]["omim_path"],
      HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
      constraint = config["trgt"]["gnomad_constraint"],
      c4r_outliers = config["trgt"]["C4R_outliers"],
      c4r = config["annotation"]["c4r"],
      promoters = config["trgt"]["promoters"],
      TR_constraint = config["trgt"]["TR_constraint"]
    resources:
        mem_mb = 20000
    conda: 
        "../envs/str_sv.yaml"
    shell: 
        """
        if [[ {params.HPO} == "none" ]]
        then
            (python3 {params.crg2_pacbio}/scripts/annotate_repeat_outliers.py --repeats {input} \
                --output_file  {output} \
                --ensembl {params.genes} \
                --gnomad_constraint {params.constraint} \
                --OMIM_path {params.OMIM} \
                --promoters {params.promoters} \
                --TR_constraint {params.TR_constraint} \
                --c4r {params.c4r} \
		        --c4r_outliers {params.c4r_outliers}) > {log} 2>&1
        else
            (python3 {params.crg2_pacbio}/scripts/annotate_repeat_outliers.py --repeats {input} \
                --output_file  {output} \
                --ensembl {params.genes} \
                --gnomad_constraint {params.constraint} \
                --OMIM_path {params.OMIM} \
                --promoters {params.promoters} \
                --TR_constraint {params.TR_constraint} \
                --c4r {params.c4r} \
                --hpo {params.HPO} \
                --c4r_outliers {params.c4r_outliers}) > {log} 2>&1
        fi
        """
          
