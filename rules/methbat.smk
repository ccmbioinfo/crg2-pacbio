rule pbcpgtools:
    input:
        samples = config["trgt"]["samples"]
    params:
        ref = config["ref"]["genome"],
        pbcpgtools = config["pbcpgtools"]["binary"], 
        model = config["pbcpgtools"]["model"] 
    output:
        combined = "pb-cpg-tools/{sample}/{sample}.combined.bed",
        hap1 = "pb-cpg-tools/{sample}/{sample}.hap1.bed",
        hap2 = "pb-cpg-tools/{sample}/{sample}.hap2.bed"
    log:
        "logs/pbcpgtools/{sample}.log"
    resources:
        mem_mb = 40000,
        threads = 8
    shell:
        """
        BAM=`cat {input.samples} | grep -P "{wildcards.sample}\t" | awk '{{print $2}}'`
        {params.pbcpgtools} \
            --threads {resources.threads} \
            --bam $BAM \
            --ref {params.ref} \
            --output-prefix pb-cpg-tools/{wildcards.sample}/{wildcards.sample} \
            --min-mapq 1 \
            --min-coverage 5 \
            --model {params.model}
        """

rule mosdepth_tiles:
    input:  
        samples = config["trgt"]["samples"]
    params:
        mosdepth = config["tools"]["mosdepth"],
        regions = config["methbat"]["regions_bed"]
    output:
        "mosdepth/{sample}/{sample}.regions.bed.gz"
    log:
        "logs/mosdepth/{sample}.log"
    resources:
        mem_mb = 40000,
        threads = 8
    shell:
        """
        BAM=`cat {input.samples} | grep -P "{wildcards.sample}\t" | awk '{{print $2}}'`
        {params.mosdepth} -t {resources.threads} -b {params.regions} -n mosdepth/{wildcards.sample}/{wildcards.sample} $BAM 
        """

rule methbat_profile:
    input: "pb-cpg-tools/{sample}/{sample}.combined.bed"
    params: 
        regions = config["methbat"]["regions"]
    output: "methbat/profiles/{sample}.profile.tsv"
    log:
        "logs/methbat/profile/{sample}.log"
    conda: "../envs/methbat.yaml"
    resources:
        mem_mb = 40000
    shell:
        """
        input=`echo {input} | sed 's/.combined.bed//'`
        methbat profile \
        --input-prefix $input \
        --input-regions {params.regions} \
        --output-region-profile {output}
        """

rule methbat_build_cohort:
    input: get_methbat_profiles,
    output: "methbat/cohort_methylation_distributions.tsv"
    log:
        "logs/methbat/build_cohort.log"
    conda: "../envs/methbat.yaml"
    resources:
        mem_mb = 600000
    shell:
        """
        # only include control samples
        echo -e "identifier\tfilename\tlabels" > methbat/cohort_profile_paths.tsv
        for sample in  {input}
        do
                ID=$(basename $sample .profile.tsv)
                echo -e "$ID\t$sample\t" >> methbat/cohort_profile_paths.tsv
        done

        methbat build \
        --input-collection  methbat/cohort_profile_paths.tsv \
        --output-profile {output}
        """


rule methbat_call_outliers:
    input:
        meth_probs = "pb-cpg-tools/{sample}/{sample}.combined.bed",
        cohort = "methbat/cohort_methylation_distributions.tsv"
    output:
        "methbat/outliers/{sample}.outliers.tsv"
    log:
        "logs/methbat/outliers/{sample}.log"
    conda: "../envs/methbat.yaml"
    resources:
        mem_mb = 40000
    shell:
        """
        input=`echo {input.meth_probs} | sed 's/.combined.bed//'`
        methbat profile \
        --input-prefix $input \
        --input-regions {input.cohort} \
        --output-region-profile {output}
        """


rule methbat_annotate_outliers:
    input: 
        outliers = "methbat/outliers/{sample}.outliers.tsv",
        coverage = "mosdepth/{sample}/{sample}.regions.bed.gz"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        ensembl = config["trgt"]["ensembl"],
        gnomad_constraint = config["trgt"]["gnomad_constraint"],
        OMIM_path = config["annotation"]["omim_path"], 
        hpo = config["run"]["hpo"]
    output: "methbat/outliers/{sample}.annotated.outliers.csv"
    log:
        "logs/methbat/annotation/{sample}.log"
    conda: "../envs/str_sv.yaml"
    resources:
        mem_mb = 40000
    shell:
        """
        if [ ! -z "{params.hpo}" ]; then
                python3 {params.crg2_pacbio}/scripts/annotate_methylation_outliers.py \
                --outliers {input.outliers} \
                --output_file {output} \
                --ensembl {params.ensembl} \
                --gnomad_constraint {params.gnomad_constraint} \
                --OMIM_path {params.OMIM_path} \
                --coverage {input.coverage} \
                --hpo {params.hpo}
        else
            echo "HPO not found"
            python3 {params.crg2_pacbio}/scripts/annotate_methylation_outliers.py \
                --outliers {input.outliers} \
                --output_file {output} \
                --ensembl {params.ensembl} \
                --gnomad_constraint {params.gnomad_constraint} \
                --OMIM_path {params.OMIM_path} \
                --coverage {input.coverage}
        fi
        """
