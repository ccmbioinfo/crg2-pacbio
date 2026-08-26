rule peddy_unphase_vcf:
    input:
        vcf=get_smallvariants_vcf
    output:
        vcf="qc/peddy/{family}.unphased.vcf.gz",
        tbi="qc/peddy/{family}.unphased.vcf.gz.tbi"
    log:
        "logs/qc/peddy/{family}.unphase.log"
    conda:
        "../envs/peddy.yaml"
    shell:
        '''
        mkdir -p qc/peddy
        bcftools view {input.vcf} \
        | sed 's/|/\//g' \
        | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

rule peddy:
    input:
        vcf=temp("qc/peddy/{family}.unphased.vcf.gz"),
        ped=format_pedigree_qc
    output:
        pca="qc/peddy/{family}.background_pca.json",
        pca_png="qc/peddy/{family}.pca_check.png",
        html="qc/peddy/{family}.html",
        het="qc/peddy/{family}.het_check.csv",
        het_png="qc/peddy/{family}.het_check.png",
        vs="qc/peddy/{family}.vs.html",
        sex="qc/peddy/{family}.sex_check.csv",
        sex_png="qc/peddy/{family}.sex_check.png",
        pedcheck="qc/peddy/{family}.ped_check.csv",
        pedcheck_png="qc/peddy/{family}.ped_check.png",
        peddy_ped="qc/peddy/{family}.peddy.ped",
        rel_diff="qc/peddy/{family}.ped_check.rel-difference.csv"
    log:
        "logs/qc/peddy/{family}.log"
    conda:
        "../envs/peddy.yaml"
    shell:
        '''
        mkdir -p qc/peddy
        peddy \
          --prefix ./qc/peddy/{wildcards.family} \
          --plot \
          --sites hg38 \
          {input.vcf} \
          {input.ped} \
          2>&1 | tee {log}
        '''

rule peddy_relatedness_mqc:
    input:
        pedcheck="qc/peddy/{family}.ped_check.csv"
    output:
        tsv="qc/multiqc_custom/{family}/peddy_relatedness_mqc.tsv"
    log:
        "logs/qc/peddy/{family}.relatedness_mqc.log"
    shell:
        '''
        mkdir -p qc/multiqc_custom/{wildcards.family}
        awk -F',' '
        BEGIN {{
            OFS="\t";
            print "Pair_ID","Sample_A","Sample_B","Peddy_Relatedness"
        }}
        NR>1 {{
            pid = $1 "_" $2
            rel = $3
            print pid, $1, $2, rel
        }}
        ' {input.pedcheck} > {output.tsv}
        '''

rule nanoplot:
    input:
        bam=get_bam
    output:
        stats="qc/nanoplot/{sample}/NanoStats.txt",
        html="qc/nanoplot/{sample}/NanoPlot-report.html",
        readlen_plot="qc/nanoplot/{sample}/WeightedHistogramReadlength.png"
    log:
        "logs/qc/nanoplot/{sample}.log"
    conda:
        "../envs/nanoplot.yaml"
    threads: 4
    shell:
        '''
        mkdir -p qc/nanoplot/{wildcards.sample}
        NanoPlot \
          --bam {input.bam} \
          -t {threads} \
          --N50 \
          --title {wildcards.sample} \
          --outdir qc/nanoplot/{wildcards.sample} \
          &> {log}
        '''

rule nanoplot_rename:
    input:
        stats="qc/nanoplot/{sample}/NanoStats.txt"
    output:
        renamed_stats="qc/nanoplot/{sample}/{sample}.txt"
    shell:
        '''
            mv {input.stats} {output.renamed_stats}
        '''

rule nanoplot_readlen:
    input:
        plot="qc/nanoplot/{sample}/WeightedHistogramReadlength.png"
    output:
        renamed_plot="qc/nanoplot/{sample}/NanoPlot_Readlength_{sample}_mqc.png"
    shell:
        '''
            mv {input.plot} {output.renamed_plot}
        '''

rule add_dp_qc:
    input:
        vcf=get_smallvariants_vcf
    output:
        temp("qc/bcftools/{family}.smallvariants_withdp.vcf")
    log:
        "logs/qc/bcftools/{family}.add_dp.log"
    wrapper:
        get_wrapper_path("bcftools", "fill-tags-qc")

rule bcftools_stats:
    input:
        vcf=f"qc/bcftools/{project}.smallvariants_withdp.vcf"
    output:
        stats="qc/bcftools/{sample}.stats"
    log: 
        "logs/qc/bcftools/{sample}.stats.log"
    conda:
        "../envs/common.yaml"
    shell:
        '''
        bcftools stats \
            -s {wildcards.sample} \
            {input.vcf} \
        | awk -v sample="{wildcards.sample}" ' 
            BEGIN {{ OFS="\t" }} 
            $1=="ID" && $2=="0" {{ $3=sample }}
            {{ print }} 
            ' \
            > {output.stats}
        '''

rule samtools_stats:
    input:
        bam=get_bam
    output:
        stats="qc/samtools/{sample}.stats"
    log:
        "logs/qc/samtools/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    params:
        ref=config["ref"]["genome"]
    shell:
        '''
        mkdir -p qc/samtools
        samtools stats -r {params.ref} {input.bam} > {output.stats} 
        '''    

rule mpile_qc:
    input:
        bam=get_bam
    output:
        pileup="qc/verifybam/{sample}.pileup"
    log:
        "logs/qc/verifybam/{sample}.mpile.log"
    conda:
        "../envs/samtools.yaml"
    params:
        ref=config["ref"]["genome"],
        svdp=config["qc"]["svdp"]
    shell:
        '''
        mkdir -p qc/verifybam
        samtools mpileup -s -B -f {params.ref} -l {params.svdp}.bed {input.bam} > {output.pileup}
        '''

rule verifybam:
    input:
        pileup="qc/verifybam/{sample}.pileup"
    output:
        selfsm="qc/verifybam/{sample}.selfSM"
    log:
        "logs/qc/verifybam/{sample}.verifybam.log"
    params:
        out_prefix="qc/verifybam/{sample}",
        sample="{sample}",
        ref=config["ref"]["genome"],
        svdp=config["qc"]["svdp"]
    wrapper:
        get_wrapper_path("verifybamid")

rule qc_pass_fail:
    input:
        selfsm=expand("qc/verifybam/{sample}.selfSM", sample=samples.index),
        sex_check="qc/peddy/{family}.sex_check.csv",
        ped_check="qc/peddy/{family}.ped_check.csv"
    output:
        tsv="qc/multiqc_custom/{family}/qc_pass_fail_mqc.tsv"
    log:
        "logs/qc/qc_pass_fail/{family}.log"
    conda:
        "../envs/str_sv.yaml"
    params:
        samples=list(samples.index),
        min_mean_coverage=config["qc"]["pass_fail_thresholds"]["min_mean_coverage"],
        max_freemix=config["qc"]["pass_fail_thresholds"]["max_freemix"],
        unrelated_max_rel=config["qc"]["pass_fail_thresholds"]["unrelated_max_rel"],
        firstdeg_min_rel=config["qc"]["pass_fail_thresholds"]["firstdeg_min_rel"],
        firstdeg_max_rel=config["qc"]["pass_fail_thresholds"]["firstdeg_max_rel"]
    script:
        "../scripts/qc_pass_fail_to_mqc.py"

rule multiqc:
    input:
        peddy_html=f"qc/peddy/{project}.html",
        peddy_relatedness="qc/multiqc_custom/{family}/peddy_relatedness_mqc.tsv",
        nanoplot_stats=expand("qc/nanoplot/{sample}/{sample}.txt", sample=samples.index),
        bcftools_stats=expand("qc/bcftools/{sample}.stats", sample=samples.index),
        selfsm=expand("qc/verifybam/{sample}.selfSM", sample=samples.index),
        samtools_stats=expand("qc/samtools/{sample}.stats", sample=samples.index),
        nanoplot_readlen=expand("qc/nanoplot/{sample}/NanoPlot_Readlength_{sample}_mqc.png", sample=samples.index),
        qc_pass_fail="qc/multiqc_custom/{family}/qc_pass_fail_mqc.tsv"
    output:
        report="qc/multiqc/{family}.multiqc_report.html"
    log:
        "logs/qc/multiqc/{family}.multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        '''
        mkdir -p qc/multiqc
        multiqc qc \
        --force \
        --filename {wildcards.family}.multiqc_report.html \
        --config {workflow.basedir}/rules/multiqc_config.yaml \
        -o qc/multiqc \
        &> {log}
        '''

rule publish_multiqc_report:
    input:
        report="qc/multiqc/{family}.multiqc_report.html"
    output:
        published="reports/{family}.multiqc_report.html"
    log:
        "logs/qc/multiqc/{family}.publish.log"
    shell:
        '''
        mkdir -p reports
        cp {input.report} {output.published}
        '''
