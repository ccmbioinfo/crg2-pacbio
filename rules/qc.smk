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
        vcf="qc/peddy/{family}.unphased.vcf.gz",
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
        stats="qc/nanoplot/{family}_{sample}/NanoStats.txt",
        html="qc/nanoplot/{family}_{sample}/NanoPlot-report.html",
        readlen_plot="qc/nanoplot/{family}_{sample}/WeightedHistogramReadlength.png"
    log:
        "logs/qc/nanoplot/{family}_{sample}.log"
    conda:
        "../envs/nanoplot.yaml"
    threads: 4
    shell:
        '''
        mkdir -p qc/nanoplot/{wildcards.family}_{wildcards.sample}
        NanoPlot \
          --bam {input.bam} \
          -t {threads} \
          --N50 \
          --title {wildcards.family}_{wildcards.sample} \
          --outdir qc/nanoplot/{wildcards.family}_{wildcards.sample} \
          &> {log}
        '''

rule nanoplot_rename:
    input:
        stats="qc/nanoplot/{family}_{sample}/NanoStats.txt"
    output:
        renamed_stats="qc/nanoplot/{family}_{sample}/{family}_{sample}.txt"
    shell:
        '''
            mv {input.stats} {output.renamed_stats}
        '''

rule nanoplot_readlen:
    input:
        plot="qc/nanoplot/{family}_{sample}/WeightedHistogramReadlength.png"
    output:
        renamed_plot="qc/nanoplot/{family}_{sample}/NanoPlot_Readlength_{family}_{sample}_mqc.png"
    shell:
        '''
            mv {input.plot} {output.renamed_plot}
        '''

rule add_dp_qc:
    input:
        vcf=get_smallvariants_vcf
    output:
        "qc/bcftools/{family}.smallvariants_withdp.vcf"
    log:
        "logs/qc/bcftools/{family}.add_dp.log"
    wrapper:
        get_wrapper_path("bcftools", "fill-tags-qc")

rule bcftools_stats:
    input:
        vcf="qc/bcftools/{family}.smallvariants_withdp.vcf"
    output:
        stats="qc/bcftools/{family}_{sample}.stats"
    log: 
        "logs/qc/bcftools/{family}_{sample}.stats.log"
    conda:
        "../envs/common.yaml"
    shell:
        '''
        bcftools stats \
            -s {wildcards.family}_{wildcards.sample} \
            {input.vcf} \
        | awk -v sample="{wildcards.family}_{wildcards.sample}" ' 
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
        stats="qc/samtools/{family}_{sample}.stats"
    log:
        "logs/qc/samtools/{family}_{sample}.log"
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
        pileup="qc/verifybam/{family}_{sample}.pileup"
    log:
        "logs/qc/verifybam/{family}_{sample}.mpile.log"
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
        pileup="qc/verifybam/{family}_{sample}.pileup"
    output:
        selfsm="qc/verifybam/{family}_{sample}.selfSM"
    log:
        "logs/qc/verifybam/{family}_{sample}.verifybam.log"
    params:
        out_prefix="qc/verifybam/{family}_{sample}",
        sample="{family}_{sample}",
        ref=config["ref"]["genome"],
        svdp=config["qc"]["svdp"]
    wrapper:
        get_wrapper_path("verifybamid")

rule multiqc:
    input:
        peddy_html=f"qc/peddy/{project}.html",
        peddy_relatedness="qc/multiqc_custom/{family}/peddy_relatedness_mqc.tsv",
        nanoplot_stats=expand("qc/nanoplot/{family}_{sample}/{family}_{sample}.txt", family=project, sample=samples.index),
        bcftools_stats=expand("qc/bcftools/{family}_{sample}.stats", family=project, sample=samples.index),
        selfsm=expand("qc/verifybam/{family}_{sample}.selfSM", family=project, sample=samples.index),
        samtools_stats=expand("qc/samtools/{family}_{sample}.stats", family=project, sample=samples.index),
        nanoplot_readlen=expand("qc/nanoplot/{family}_{sample}/NanoPlot_Readlength_{family}_{sample}_mqc.png", family=project, sample=samples.index)
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