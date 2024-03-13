rule generate_allele_db:
    input: 
        vcf_path = config["trgt"]["vcf_path"]
    output: 
        out_path = "allele_db/alleles_db.gz"
    log: "logs/allele_db/allele_db.log"
    resources:
        mem_mb = 10000
    conda:
        "../envs/outlier_expansions.yaml"
    script:
        "../scripts/generate_allele_db.py"

rule sort_allele_db:
    input: 
        input_dir = "allele_db/",
        input_file = "allele_db/alleles_db.gz"
    output: "allele_db/sorted_alleles_db.gz"
    log: "logs/allele_db/sort_allele_db.log"
    shell:
        """
        (zcat < {input.input_file} | sort -T {input.input_dir}  -k 1,1 | gzip > {output}) > {log} 2>&1
        """

rule find_repeat_outliers:
    input: 
        alleles_path = "allele_db/sorted_alleles_db.gz",
        case_ids = config["trgt"]["samples"]
    output:
        out_path =  "repeat_outliers/repeat_outliers.csv"
    log:  "logs/repeat_outliers/repeat_outliers.log"
    resources:
        mem_mb = 20000
    conda: 
        "../envs/outlier_expansions.yaml"
    script:
        "../scripts/find_repeat_outliers.py"


rule annotate_repeat_outliers:
    input: "repeat_outliers/repeat_outliers.csv"
    output: "repeat_outliers/annotated_repeat_outliers.csv"
    log:  "logs/repeat_outliers/annotate_repeat_outliers.log"
    params: 
      crg2_path = config["run"]["crg2_path"],
      genes = config["trgt"]["ensembl_gtf"],
      OMIM = config["trgt"]["omim_path"],
      HPO = config["trgt"]["hpo"],
      constraint = config["trgt"]["gnomad_constraint"]
    resources:
        mem_mb = 20000
    conda: 
        "../envs/outlier_expansions.yaml"
    shell: 
        """
        (python3 {params.crg2_path}/scripts/annotate_repeat_outliers.py --repeats {input} \
            --output_file  {output} \
            --ensembl_gtf {params.genes} \
            --gnomad_constraint {params.constraint} \
            --OMIM_path {params.OMIM} \
            --hpo {params.HPO}) > {log} 2>&1
        """
          
