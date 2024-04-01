rule generate_allele_db:
    input: config["run"]["units"]
    output: "repeat_outliers/alleles_db.gz"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        family = config["run"]["project"]
    log: "logs/repeat_outliers/allele_db.log"
    resources:
        mem_mb = 10000
    conda:
        "../envs/outlier_expansions.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/generate_allele_db.py --vcf_path {input} \
            --output_file  {output} --family {params.family}) > {log} 2>&1
        """

rule sort_allele_db:
    input: 
        input_file = "repeat_outliers/alleles_db.gz"
    output: "repeat_outliers/sorted_alleles_db.gz"
    log: "logs/repeat_outliers/sort_allele_db.log"
    shell:
        """
        (zcat < {input.input_file} | sort -T . -k 1,1 | gzip > {output}) > {log} 2>&1
        """

rule find_repeat_outliers:
    input: 
        alleles_path = "repeat_outliers/sorted_alleles_db.gz",
        case_ids = config["trgt"]["samples"]
    output:
        out_path =  "repeat_outliers/repeat_outliers.csv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    log:  "logs/repeat_outliers/repeat_outliers.log"
    resources:
        mem_mb = 20000
    conda: 
        "../envs/outlier_expansions.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/find_repeat_outliers.py --alleles_path {input.alleles_path} \
            --output_file  {output} \
            --case_ids {input.case_ids}) > {log} 2>&1
        """


rule annotate_repeat_outliers:
    input: "repeat_outliers/repeat_outliers.csv"
    output: "repeat_outliers/annotated_repeat_outliers.csv"
    log:  "logs/repeat_outliers/annotate_repeat_outliers.log"
    params: 
      crg2_pacbio = config["tools"]["crg2_pacbio"],
      genes = config["trgt"]["ensembl_gtf"],
      OMIM = config["trgt"]["omim_path"],
      HPO = config["trgt"]["hpo"] if config["trgt"]["hpo"] else "none",
      constraint = config["trgt"]["gnomad_constraint"]
    resources:
        mem_mb = 20000
    conda: 
        "../envs/outlier_expansions.yaml"
    shell: 
        """
        if [[ {params.HPO} == "none" ]]
        then
            (python3 {params.crg2_pacbio}/scripts/annotate_repeat_outliers.py --repeats {input} \
                --output_file  {output} \
                --ensembl_gtf {params.genes} \
                --gnomad_constraint {params.constraint} \
                --OMIM_path {params.OMIM}) > {log} 2>&1
        else
            (python3 {params.crg2_pacbio}/scripts/annotate_repeat_outliers.py --repeats {input} \
                --output_file  {output} \
                --ensembl_gtf {params.genes} \
                --gnomad_constraint {params.constraint} \
                --OMIM_path {params.OMIM} \
                --hpo {params.HPO}) > {log} 2>&1
        fi
        """
          
