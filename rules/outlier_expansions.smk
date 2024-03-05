rule generate_allele_db:
    input = config["run"]["TRGT_VCF_path"],
    output =  "allele_db/alleles_db.gz",
    log = "logs/allele_db/allele_db.log"
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "../envs/trgt.yaml"
    script:
        "../scripts/generate_allele_db.py"

rule sort_allele_db:
    input:
    output:
    log:
    params:
    resources:
    conda: 
    shell:

rule extract_repeat_alleles:
    input:
    output:
    log:
    params:
    resources:
    conda: 
    script:

rule find_repeat_outliers:
    input:
    output:
    log:
    params:
    resources:
    conda: 
    script:

rule annotate_repeat_outliers:
    input:
    output:
    log:
    params:
    resources:
    conda: 
    script: