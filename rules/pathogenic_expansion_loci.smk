# rule genotype_pathogenic_loci:
# /hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread/TRGT/scripts/genotype_repeats.sh

rule merge_vcfs:
    input: get_trgt_path_str_vcf_dir
    output: "pathogenic_repeats/{family}.known.path.str.loci.vcf"
    log: "logs/pathogenic_repeats/{family}.merge.repeat.vcfs.log"
    conda:
        "../envs/common.yaml"
    shell: 
        """
        echo {input}/*vcf.gz | tr ' ' '\\n' > pathogenic_repeats/{wildcards.family}_vcfs.txt
        bcftools merge -l pathogenic_repeats/{wildcards.family}_vcfs.txt -o pathogenic_repeats/{wildcards.family}.known.path.str.loci.vcf -O v -F x -m all
        """

rule annotate_pathogenic_repeats:
    input: "pathogenic_repeats/{family}.known.path.str.loci.vcf"
    output: "pathogenic_repeats/{family}.known.path.str.loci.csv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        disease_thresholds = config["annotation"]["pathogenic_repeats"]["disease_thresholds"]
    log: "logs/pathogenic_repeats/{family}.annotate.repeats.log"
    resources:
        mem_mb = 10000
    conda:
        "../envs/str_sv.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/annotate_path_str_loci.py --vcf {input} \
            --output_file  {output} \
            --disease_thresholds {params.disease_thresholds}) > {log} 2>&1
        """