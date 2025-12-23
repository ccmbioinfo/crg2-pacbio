rule pbcpgtools:
    input:
        samples = config["run"]["samples"]
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
        samples = config["run"]["samples"]
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

rule filter_sample_smvs: 
    input: get_sample_smvs
    output: 
        temp("methbat/{sample}.rare.smvs.bed")
    conda: "../envs/common.yaml"
    log:
        "logs/methbat/filter_sample_smvs/{sample}.log"
    shell:
        """
        {{
            bcftools query -s {wildcards.sample} -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
        }} || {{
            # renamed genesteps samples 
            family=`echo {wildcards.sample}| cut -d'_' -f1,2,3 | tr -d '_'`
            participant=`echo {wildcards.sample}| cut -d'_' -f4`
            sample=`echo -e "${{family}}_${{participant}}"`
            bcftools query -s $sample -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
        }} || {{
           # some samples from TCAG have RLGS suffix
           sample="{wildcards.sample}_RLGS"
           bcftools query -s $sample -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
        }} || {{
	       sample="{wildcards.sample}_RLGA"
           bcftools query -s $sample -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
        }} || {{
           participant=`echo {wildcards.sample} | cut -d'_' -f2`
           sample=`ls /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/data/ | grep -v RNA | grep $participant`
           bcftools query -s $sample -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
        }} || {{
           sample=`echo {wildcards.sample} | sed 's/_DNA_//' | sed 's/_A1//'`
           bcftools query -s $sample -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
}} || {{
           sample=`echo {wildcards.sample} | sed 's/DSK_/DSK/'`
           bcftools query -s $sample -i '(INFO/gnomad_af_popmax <0.01 | INFO/gnomad_af_popmax == ".") & GT!="RR" & GT!="./." ' {input} -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t[%GT]\t%gnomad_af_popmax\n' > {output}
}}
        """

rule methbat_annotate_outliers:
    input: 
        outliers = "methbat/outliers/{sample}.outliers.tsv",
        coverage = "mosdepth/{sample}/{sample}.regions.bed.gz",
        smvs ="methbat/{sample}.rare.smvs.bed",
        svs = get_pbsv_csv,
        cnvs = get_cnvs,
        trs = get_TR_outliers
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        ensembl = config["annotation"]["general"]["ensembl"],
        gnomad_constraint = config["annotation"]["general"]["gnomad_constraint"],
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
                --SV_path {input.svs} \
                --CNV_path {input.cnvs} \
                --TR_outlier_path {input.trs} \
                --SMV_path {input.smvs} \
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
                --SV_path {input.svs} \
                --CNV_path {input.cnvs} \
                --TR_outlier_path {input.trs} \
                --SMV_path {input.smvs} \
                --coverage {input.coverage} 
        fi
        """
