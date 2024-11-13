def get_filt_vcf(wildcards):
    if wildcards.p == "coding":
        return "filtered/{family}.vcf.gz"
    elif wildcards.p == "denovo":
        return "filtered/{family}.vcf.gz"
    else:
        return "filtered/{p}/{family}.{p}.vcf.gz".format(p=wildcards.p,family=project)


rule input_prep:
    input:
        units=config["run"]["units"]
    params:
        outdir= "filtered"
    output:
        "filtered/{family}.vcf.gz"
    wildcard_constraints:
        family = "(?!.*panel|.*coding).*"
    log:
        "logs/input_prep/{family}.log"
    threads:
        4
    wrapper:
        get_wrapper_path("input_prep")

rule vt:
    input: get_filt_vcf # (vcf, bcf, or vcf.gz)
    output:
        temp("filtered/{p}/{family}.{p}.uniq.normalized.decomposed.vcf"),
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{family}.vt.{p}.uniq.normalized.decomposed.log"
    wrapper:
        get_wrapper_path("vt")

rule pass:
    input:
       	"{prefix}.{ext}"
    output:
        temp("{prefix}.pass.{ext,(vcf|vcf\.gz)}")
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: "-f PASS"
    wrapper:
        get_wrapper_path("bcftools", "view")

rule vep:
    input:
        "filtered/{p}/{family}.{p}.uniq.normalized.decomposed.pass.vcf",
    output:
        temp("annotated/{p}/vep/{family}.{p}.vep.vcf"),
    log:
        "logs/vep/{family}.vep.{p}.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        ref=config["ref"]["genome"],
    wrapper:
        get_wrapper_path("vep")

rule vcfanno:
    input:
        "annotated/{p}/vep/{family}.{p}.vep.vcf",
    output:
        temp("annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.noDP.noPS_info.vcf.gz"),
    log:
        "logs/vcfanno/{family}.vcfanno.{p}.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
        lua_script=config["annotation"]["pacbio.vcfanno"]["lua_script"],
       	conf=config["annotation"]["pacbio.vcfanno"]["conf"],
        base_path=config["annotation"]["pacbio.vcfanno"]["base_path"],
    wrapper:
        get_wrapper_path("vcfanno")

rule add_dp_field:
    input: 
        "annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.noDP.noPS_info.vcf.gz",
    output:
        temp("annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.noPS_info.vcf"),
    log:
        "logs/bcftools/{family}.add_dp_field.{p}.log"
    wrapper:
        get_wrapper_path("bcftools","fill-tags")

rule add_ps_field:
    input:
       vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.noPS_info.vcf",
    output:
        temp("annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf"),
    log:
        "logs/bcftools/{family}.add_PS_INFO_field.{p}.log"
    shell:
        '''
            #Get the PS values for every sample from the FORMAT VCF field, remove the trailing comma and store the results in PS_annot.txt
            bcftools query -f '%CHROM\t%POS\t[%PS,]\n' {input.vcf} | sed 's/,*$//g' > annotated/{wildcards.p}/vcfanno/PS_annot.txt

            #BGZIP and TABIX the PS_annot.txt file
            bgzip annotated/{wildcards.p}/vcfanno/PS_annot.txt
            tabix -s1 -b2 -e2 annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz

            #create header file containing PS INFO field info
            echo -e "##INFO=<ID=PS,Number=.,Type=String,Description="Phase set">" > annotated/{wildcards.p}/vcfanno/hdr.txt

            #Annotate the VCF with the PS_annot.txt.gz file to add PS tag info to the VCF
            bcftools annotate -a annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz -h annotated/{wildcards.p}/vcfanno/hdr.txt -c CHROM,POS,INFO/PS {input.vcf} > {output}

            #Remove intermediate files  
            rm annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz.tbi annotated/{wildcards.p}/vcfanno/hdr.txt
        '''

rule vcf2db:
    input:
        "annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf",
    output:
         db="annotated/{p}/{family}-gemini.db",
    log:
        "logs/vcf2db/{family}.vcf2db.{p}.log"
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")


rule bgzip:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    wildcard_constraints:
        prefix = "(?!.*trgt).*"
    conda:
        "../envs/common.yaml"

    shell:
        '''
        bgzip -c {input} > {output}
        '''

rule tabix:
    input: 
        "{prefix}.vcf.gz"
    output: 
        "{prefix}.vcf.gz.tbi"
    log: 
        "logs/{prefix}.log"
    conda:
           "../envs/common.yaml"
    wrapper:
        get_wrapper_path("tabix")
