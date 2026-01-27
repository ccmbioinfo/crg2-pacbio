rule genotype_adotto_loci:
    input: get_bam
    output: 
        vcf = temp("repeat_outliers/{family}_{sample}.trgt.unsorted.vcf.gz"),
        spanning_bam = "repeat_outliers/{family}_{sample}.trgt.unsorted.spanning.bam"
    params: 
        trgt = config["tools"]["trgt"],
        ref = config["ref"]["genome"],
        repeats = config["annotation"]["general"]["adotto_repeats"]
    log: "logs/repeat_outliers/{family}_{sample}.trgt.log"
    conda:
        "../envs/common.yaml"
    threads: 10
    resources:
        mem_mb = 60000
    shell: 
        """
        # get sex, code from https://github.com/ccmbioinfo/crg/get_XY.sh
	module load samtools
        x=`samtools idxstats {input} | egrep  "X|chrX"`;
        y=`samtools idxstats {input} | egrep  "Y|chrY"`;
        xcov=`echo $x | awk '{{ printf("%0.5f", $3/$2); }}'`;
        ycov=`echo $y | awk '{{ printf("%0.5f", $3/$2); }}'`;

        rat=$(echo "scale=4; ${{xcov}}/${{ycov}}" | bc)
        if (( $(echo "$rat > 5.0" | bc -l) )); then
            sex=XX
        else
            sex=XY
        fi

        echo "$sex"

        {params.trgt} genotype --genome {params.ref} \
            --reads {input} \
            --repeats {params.repeats} \
            --output-prefix repeat_outliers/{wildcards.family}_{wildcards.sample}.trgt.unsorted \
            --karyotype $sex \
            --sample-name {wildcards.family}_{wildcards.sample} \
            --threads {threads}
        """       

rule sort_trgt_adotto_vcf:
    input: "repeat_outliers/{family}_{sample}.trgt.unsorted.vcf.gz"
    output: "repeat_outliers/{family}_{sample}.trgt.sorted.vcf.gz"
    log: "logs/bcftools/{family}_{sample}.sort.trgt.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        bcftools sort -Ob -o {output} {input};
        bcftools index {output}
        """

rule merge_trgt_vcf:
    input: 
        vcf = expand("repeat_outliers/{{family}}_{sample}.trgt.sorted.vcf.gz", sample=samples.index),
        indices = expand("repeat_outliers/{{family}}_{sample}.trgt.sorted.vcf.gz.tbi", sample=samples.index) 
    params:
        trgt = config["tools"]["trgt"],
        genome = config["ref"]["genome"]
    output: "repeat_outliers/{family}.trgt.vcf.gz"
    log: "logs/repeat_outliers/{family}.merge.trgt.vcf.log"
    shell:
        """ 
        export TMPDIR=`pwd`
        {params.trgt} -vv merge \
            --vcf {input.vcf} \
            --force-single \
            --genome {params.genome} \
            --output-type z > {output}
        """

rule calculate_lps: 
    input: "repeat_outliers/{family}_{sample}.trgt.sorted.vcf.gz"
    output: temp("repeat_outliers/{family}_{sample}.trgt.lps.tsv")
    params:
        trgt_lps = config["tools"]["trgt-lps"]
    log: "logs/repeat_outliers/{family}_{sample}.trgt.lps.log"
    shell:
        """
        {params.trgt_lps} --vcf {input} > {output}
        """

rule combine_lps: 
    input: 
        lps = expand("repeat_outliers/{{family}}_{sample}.trgt.lps.tsv", sample=samples.index)
    output: "repeat_outliers/{family}.trgt.lps.combined.tsv.gz"
    log: "logs/repeat_outliers/{family}.trgt.lps.combined.log"
    run:
        import pandas as pd
        files = input.lps
        lps_list = []
        for file in files: 
            print(files)
            df = pd.read_csv(file, sep="\t")
            sample = file.split("/")[-1].split(".")[0]
            df["sample"] = sample
            lps_list.append(df)
        lps_df = pd.concat(lps_list)
        lps_df.to_csv(output[0], sep="\t", index=False, compression="gzip")

# rule sort_merged_trgt_vcf:
#     input: "repeat_outliers/{family}.trgt.vcf.gz"
#     output: "repeat_outliers/{family}.trgt.sorted.vcf.gz"
#     log: "logs/repeat_outliers/{family}.trgt.sorted.vcf.log"
#     conda:
#         "../envs/common.yaml"
#     shell:
#         """
#         bcftools sort -O z -o {output} {input}
#         tabix {output}
#         """

rule find_repeat_outliers:
    input: 
        case_vcf = "repeat_outliers/{family}.trgt.vcf.gz",
        control_vcf = config["trgt"]["control_alleles"], 
        lps = "repeat_outliers/{family}.trgt.lps.combined.tsv.gz",
        control_lps = config["trgt"]["control_lps"]
    output: "reports/{family}.repeat.outliers.tsv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    log:  "logs/repeat_outliers/{family}.repeat.outliers.log"
    resources:
        mem_mb = 300000
    conda: 
        "../envs/str_sv.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/find_repeat_outliers.py --case_vcf {input.case_vcf} \
            --control_vcf {input.control_vcf} \
            --case_lps {input.lps} \
            --control_lps {input.control_lps} \
            --output_file  {output}) > {log} 2>&1 
        """

rule annotate_repeat_outliers:
    input: "repeat_outliers/{family}.repeat.outliers.tsv"
    output: "repeat_outliers/{family}.repeat.outliers.annotated.csv"
    log:  "logs/repeat_outliers/{family}.annotate.repeat.outliers.log"
    params: 
      crg2_pacbio = config["tools"]["crg2_pacbio"],
      genes = config["annotation"]["general"]["ensembl"],
      OMIM = config["annotation"]["omim_path"],
      HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
      constraint = config["annotation"]["general"]["gnomad_constraint"],
      c4r_outliers = config["trgt"]["C4R_outliers"],
      c4r = config["annotation"]["c4r"],
      promoters = config["annotation"]["general"]["promoters"],
      TR_constraint = config["trgt"]["TR_constraint"],
      repeat_catalog = config["annotation"]["general"]["adotto_repeats"]
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
                --repeat_catalog {params.repeat_catalog} \
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
                --repeat_catalog {params.repeat_catalog} \
                --hpo {params.HPO} \
                --c4r_outliers {params.c4r_outliers}) > {log} 2>&1
        fi
        """
          
