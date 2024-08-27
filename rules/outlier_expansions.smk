rule generate_allele_db:
    input: config["run"]["units"]
    output: "repeat_outliers/allele_db/{family}.alleles.db.gz"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        family = config["run"]["project"]
    log: "logs/repeat_outliers/{family}_allele_db.log"
    resources:
        mem_mb = 10000
    conda:
        "../envs/str_sv.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/generate_allele_db.py --vcf_path {input} \
            --output_file  {output} --family {params.family}) > {log} 2>&1
        """

rule sort_allele_db:
    input: 
        input_file = "repeat_outliers/allele_db/{family}.alleles.db.gz"
    output: "repeat_outliers/allele_db/{family}.alleles.sorted.db.gz"
    log: "logs/repeat_outliers/{family}.allele.sorted.db.log"
    shell:
        """
        (zcat < {input.input_file} | sort -T . -k 1,1 | gzip > {output}) > {log} 2>&1
        """

rule find_repeat_outliers:
    input: 
        alleles_path = "repeat_outliers/allele_db/{family}.alleles.sorted.db.gz",
        control_alleles = config["trgt"]["control_alleles"]
    output: "repeat_outliers/{family}.repeat.outliers.csv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    log:  "logs/repeat_outliers/{family}.repeat.outliers.log"
    resources:
        mem_mb = 100000
    conda: 
        "../envs/str_sv.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/find_repeat_outliers.py --alleles_path {input.alleles_path} \
            --control_alleles {input.control_alleles} \
            --output_file  {output}) > {log} 2>&1 
        """


rule annotate_repeat_outliers:
    input: "repeat_outliers/{family}.repeat.outliers.csv"
    output: "repeat_outliers/{family}.repeat.outliers.annotated.csv"
    log:  "logs/repeat_outliers/{family}.annotate.repeat.outliers.log"
    params: 
      crg2_pacbio = config["tools"]["crg2_pacbio"],
      genes = config["trgt"]["ensembl_gtf"],
      OMIM = config["trgt"]["omim_path"],
      HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
      constraint = config["trgt"]["gnomad_constraint"]
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
          
