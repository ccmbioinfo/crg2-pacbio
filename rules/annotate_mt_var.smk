rule mitorsaw_call:
    input:
        bam=get_bam
    output:
        vcf="mitochondrial_variants/{family}/{sample}.mitorsaw.vcf.gz",
        hap_stats="mitochondrial_variants/{family}/{sample}.hap_stats.json"
    params:
        ref=config["ref"]["genome"]
    log:
        "logs/mito/mitorsaw/{family}_{sample}.log"
    conda:
        "../envs/mitorsaw.yaml"
    shell:
        """
        mkdir -p mitochondrial_variants/{wildcards.family}/{wildcards.sample}/debug
        mitorsaw haplotype \
            --reference {params.ref} \
            --bam {input.bam} \
            --output-vcf {output.vcf} \
            --output-hap-stats {output.hap_stats} \
            --output-debug mitochondrial_variants/{wildcards.family}/{wildcards.sample}/debug \
            > {log} 2>&1
        """


rule merge_mitorsaw_vcfs:
    input:
        vcfs=expand("mitochondrial_variants/{{family}}/{sample}.mitorsaw.vcf.gz", sample=samples.index),
        indices=expand("mitochondrial_variants/{{family}}/{sample}.mitorsaw.vcf.gz.tbi", sample=samples.index)
    output:
        "mitochondrial_variants/{family}/{family}.mitorsaw.vcf.gz"
    log:
        "logs/mito/merge/{family}.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        (
        n_vcfs=$(printf '%s\n' {input.vcfs} | wc -l)
        if [ "$n_vcfs" -eq 1 ]; then
            cp {input.vcfs} {output}
            tabix -f {output}
        else
            bcftools merge -m none {input.vcfs} -Oz -o {output}
            tabix -f {output}
        fi
        ) > {log} 2>&1
        """


rule mitorsaw_normalize:
    input:
        "mitochondrial_variants/{family}/{family}.mitorsaw.vcf.gz"
    output:
        temp("mitochondrial_variants/{family}/{family}.mt.normalise.decompose.vcf.gz")
    params:
        outdir="mitochondrial_variants/{family}",
        tool=config["tools"]["mity"],
        crg2_pacbio=config["tools"]["crg2_pacbio"]
    log:
        "logs/mito/normalize/{family}.log"
    wrapper:
        get_wrapper_path("mito", "normalize")


rule mity_report:
    input:
        "mitochondrial_variants/{family}/{family}.mt.normalise.decompose.vcf.gz"
    output:
        "mitochondrial_variants/{family}/{family}.mity.annotated.vcf.gz",
        "mitochondrial_variants/{family}/{family}.mity.report.xlsx"
    params:
        outdir="mitochondrial_variants/{family}",
        tool=config["tools"]["mity"],
        report_config=config["annotation"]["mity"]["report_config"],
        vcfanno_config=config["annotation"]["mity"]["vcfanno_config"],
        base_path=config["annotation"]["mity"]["base_path"]
    log:
        "logs/mito/mity_report/{family}.log"
    wrapper:
        get_wrapper_path("mito", "report")


rule generate_mt_report:
    input:
        vcf="mitochondrial_variants/{family}/{family}.mity.annotated.vcf.gz",
        report="mitochondrial_variants/{family}/{family}.mity.report.xlsx"
    output:
        "reports/{family}.mito.csv"
    log:
        "logs/report/mitochondrial/{family}.mitochondrial.report.log"
    conda:
        "../envs/mt_report.yaml"
    script:
        "../scripts/mt_report.py"
