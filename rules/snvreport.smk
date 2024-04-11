rule allsnvreport:
    input:
        db="annotated/{p}/{family}-gemini.db",
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz"
    output:
        directory("small_variants/{p}/{family}")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/{p}/{family}.cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         database_path=config["annotation"]["cre"]["database_path"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p {output}
         cd {output}
         ln -s ../../../{input.db} {project}-ensemble.db
         #bgzip ../../../{input.vcf} -c > {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s ../../../{input.vcf} {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}-ensemble-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {project}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         if [ {wildcards.p} == "coding" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} {params.cre}/cre.sh {project} 
         elif [ {wildcards.p} == "denovo" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} type=denovo {params.cre}/cre.sh {project} 
         else
         cre={params.cre} reference={params.ref} database={params.database_path} type=wgs {params.cre}/cre.sh {project}
         unset type
         fi;
         '''