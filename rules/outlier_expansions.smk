rule genotype_adotto_loci:
    input: get_bam
    output: 
        vcf = temp("repeat_outliers/{sample}.trgt.unsorted.vcf.gz"),
        spanning_bam = temp("repeat_outliers/{sample}.trgt.unsorted.spanning.bam")
    params: 
        trgt = config["tools"]["trgt"],
        ref = config["ref"]["genome"],
        repeats = config["annotation"]["general"]["TRExplorer_catalog"]
    log: "logs/repeat_outliers/{sample}.trgt.log"
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
            --output-prefix repeat_outliers/{wildcards.sample}.trgt.unsorted \
            --karyotype $sex \
            --sample-name {wildcards.sample} \
            --threads {threads}
        """       

rule sort_trgt_adotto_vcf:
    input: "repeat_outliers/{sample}.trgt.unsorted.vcf.gz"
    output: temp("repeat_outliers/{sample}.trgt.sorted.vcf.gz")
    log: "logs/bcftools/{sample}.sort.trgt.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        bcftools sort -Ob -o {output} {input};
        bcftools index {output}
        """

rule merge_trgt_vcf:
    input: 
        vcf = expand("repeat_outliers/{sample}.trgt.sorted.vcf.gz", sample=samples.index),
        indices = expand("repeat_outliers/{sample}.trgt.sorted.vcf.gz.tbi", sample=samples.index)
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

rule filter_trgt_vcf:
    input: "repeat_outliers/{family}.trgt.vcf.gz"
    output: "repeat_outliers/{family}.trgt.filtered.vcf.gz"
    log: "logs/repeat_outliers/{family}.trgt.filtered.vcf.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        # filter out sites with no ALT allele
        bcftools view -e 'ALT="."' -Oz -o {output} {input}
        """

rule calculate_lps: 
    input: 
        vcf = "repeat_outliers/{family}.trgt.filtered.vcf.gz",
        index = "repeat_outliers/{family}.trgt.filtered.vcf.gz.tbi"
    output: temp("repeat_outliers/{family}.trgt.lps.tsv")
    params:
        trgt_lps = config["tools"]["trgt-lps"],
        threads = 8
    log: "logs/repeat_outliers/{family}.trgt.lps.log"
    shell:
        """
        {params.trgt_lps} --vcf {input.vcf} --threads {params.threads} > {output}
        """

family = config["run"]["project"]
hpo_matches = f"hpo/{family}.matched.tsv"
hpo_available = config["run"].get("hpo", "")
ped_available = config["run"].get("ped", "")

def get_annotate_repeat_outliers_args(wildcards, input):
    args = []
    if hpo_available:
        args.append(f'--hpo "{hpo_matches}"')
    if ped_available:
        args.append(f'--ped "{ped_available}"')

    if len(args) > 0:
        return ' '.join(args)
    else:
        return ""

rule annotate_repeat_outliers:
    input:
      TRGT_vcf = "repeat_outliers/{family}.trgt.filtered.vcf.gz",
      LPS = "repeat_outliers/{family}.trgt.lps.tsv",
      hpo = [hpo_matches] if config["run"]["hpo"] else [],
    output: "reports/{family}.repeat.outliers.annotated.csv"
    log:  "logs/repeat_outliers/{family}.annotate.repeat.outliers.log"
    params: 
      HPRC_AL_DBSCAN_cutoffs = config["trgt"]["HPRC_AL_DBSCAN_cutoffs"],
      HPRC_LPS_DBSCAN_cutoffs = config["trgt"]["HPRC_LPS_DBSCAN_cutoffs"],
      HPRC_stats = config["trgt"]["HPRC_stats"],
      TRExplorer_annotated_catalog = config["trgt"]["TRExplorer_annotated_catalog"],
      crg2_pacbio = config["tools"]["crg2_pacbio"],
      gencode_gtf = config["annotation"]["general"]["gencode_gtf"],
      OMIM = config["annotation"]["omim_path"],
      constraint = config["annotation"]["general"]["gnomad_constraint"],
      promoters = config["annotation"]["general"]["promoters"],
      promoter_size = config["annotation"]["general"]["promoter_size"],
      hpo_ped = get_annotate_repeat_outliers_args
    resources:
        mem_mb = 20000
    conda: 
        "../envs/str_sv.yaml"
    shell: 
        """
            (python3 {params.crg2_pacbio}/scripts/annotate_repeat_outliers.py --vcf {input.TRGT_vcf} \
                --lps {input.LPS} \
                --dbscan_cutoffs_lps {params.HPRC_LPS_DBSCAN_cutoffs} \
                --dbscan_cutoffs_al {params.HPRC_AL_DBSCAN_cutoffs} \
                --hprc_stats {params.HPRC_stats} \
                --annotated_catalog {params.TRExplorer_annotated_catalog} \
                --gtf {params.gencode_gtf} \
                --SCREEN_promoters {params.promoters} \
                --promoter_size {params.promoter_size} \
                --omim {params.OMIM} \
                {params.hpo_ped} \
                --outfile  {output}) > {log} 2>&1

        """
          
