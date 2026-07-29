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
         set -eo pipefail
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
         elif [ {wildcards.p} == "wgs-high-impact" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} type=wgs.high.impact {params.cre}/cre.sh {project}
         else
         cre={params.cre} reference={params.ref} database={params.database_path} type=wgs {params.cre}/cre.sh {project}
         unset type
         fi;
         '''

if config["run"]["hpo"]:

    def get_panel(wildcards):
        if wildcards.p == "panel":
            return "genes/{family}.bed"
        else:
            return "genes/{family}_{p}.bed"
        return config["run"]["panel"]

    # def get_bed(wildcards):
    #     if wildcards.p == "panel-flank":
    #         return "genes/{family}_{p}.bed"
    #     return get_panel()


    rule hpo_to_panel:
        input: 
            hpo=config["run"]["hpo"],
            ensembl=config["genes"]["ensembl"],
            refseq=config["genes"]["refseq"],
            hgnc=config["genes"]["hgnc"]
        params: 
            crg2_pacbio=config["tools"]["crg2_pacbio"],
            cre=config["tools"]["cre"]
        output: 
            genes="genes/{family}.bed"
        wildcard_constraints:
            family = "(?!.*panel|.*coding).*"
        conda: "../envs/hpo_to_panel.yaml"
        log: "logs/hpo_to_panel/{family}.log"
        script:
            "../scripts/hpo_to_panel.py"

    rule add_flank:
        input: "genes/{family}.bed"
        output: "genes/{family}_{p}.bed"
        params: config["run"]["flank"]
        shell:
            '''
            cat {input} | awk -F "\t" '{{print $1"\t"$2-{params}"\t"$3+{params}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output}
            '''

    rule intersect:
        input: 
            left="filtered/{family}.vcf.gz",
            right=get_panel
        output:
            vcf="filtered/{p}/{family}.{p}.vcf.gz"
        params:
            extra="-header"
        log: "logs/report/bedtools-{family}-{p}.log"
        wrapper:
            get_wrapper_path("bedtools", "intersect")
            
        
    rule annotate_hpo:
        input:
            reports=expand("small_variants/{p}/{family}",p=["coding", "panel", "panel-flank", "denovo"] if config["run"]["ped"] else ["coding", "panel", "panel-flank"], family=project),
            hpo=config["run"]["hpo"]
        output: 
            directory("report/hpo_annotated")
        conda: 
            "../envs/hpo_to_panel.yaml"
        params: 
            crg2_pacbio = config["tools"]["crg2_pacbio"]
        log: 
            "logs/hpo_annotation.log"
        shell:
            '''
                if [ ! -d {output} ]; then mkdir -p {output}; fi;
                for i in {input.reports}; do 
                    echo "dir: ${{i}}" >> {log};                
                    if [[ "${{i}}" =~ .*"panel-flank".* ]]; then 
                        j=`find ${{i}} -name  "*.wgs.[0-9]*.csv" | grep -v "clinical"`;
                        rename=`echo ${{j}} | sed 's/wgs/wgs.panel-flank100k/g'`;
                        if [ ! -f ${{rename}} ]; then
                            ln -s `basename ${{j}}` ${{rename}};
                        else
                            echo "${{j}} found panel-flank"  >> {log};
                        fi;
                    elif [[ "${{i}}" =~ .*"panel".* ]]; then 
                        j=`find ${{i}} -name "*.wgs.[0-9]*.csv" | grep -v "clinical"`;
                        rename=`echo ${{j}} | sed 's/wgs/wgs.panel/g'`;
                        if [ ! -f ${{rename}} ]; then
                            ln -s `basename ${{j}}` ${{rename}};
                        else
                            echo "${{j}} found panel"  >> {log};
                        fi;
                    elif [[ "${{i}}" =~ .*"denovo".* ]]; then 
                        rename=`find ${{i}} -name  "*.denovo.[0-9]*.csv" | grep -v "clinical"`;
                        echo "${{rename}} found denovo"  >> {log};
                    else 
                        rename=`find ${{i}} -name "*.wes*.[0-9]*.csv" | grep -v "clinical"`;
                        echo "${{rename}} found wes"  >> {log};
                    fi;
                    python {params.crg2_pacbio}/scripts/add_hpo_terms_to_wes.py {input.hpo} ${{rename}} {output} >> {log} 2>&1
                done;
            '''


rule slivar_select:
    input:
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz"
    output:
        rare_main=temp("small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz"),
        rare_main_tbi=temp("small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz.tbi"),
        rare_clinvar=temp("small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz"),
        rare_clinvar_tbi=temp("small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz.tbi"),
        common_pathogenic_clinvar=temp("small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz"),
        common_pathogenic_clinvar_tbi=temp("small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz.tbi"),
    wildcard_constraints:
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.select.log"
    conda:
        os.path.join(os.path.expanduser(config["tools"]["cphi_dragen_anno"]), "workflow", "envs", "slivar.yaml")
    params:
        js=config["tools"]["cphi_dragen_anno"] + "/workflow/scripts/slivar/slivar_functions.js",
        consequence_order_file=config["tools"]["cphi_dragen_anno"] + "/workflow/scripts/slivar/default-order.txt",
        mode="{p}",
        profile="pacbio",
    wrapper:
        "file:" + os.path.join(os.path.expanduser(config["tools"]["cphi_dragen_anno"]), "workflow", "wrappers", "slivar")


rule slivar_postfilter:
    input:
        rare_main="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz",
        rare_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz",
        common_pathogenic_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz",
    output:
        vcf=temp("small_variants_slivar/{p}/{family}/{family}.{p}.postfilter.vcf"),
    wildcard_constraints:
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.postfilter.log"
    conda:
        os.path.join(os.path.expanduser(config["tools"]["cphi_dragen_anno"]), "workflow", "envs", "slivar.yaml")
    params:
        cphi_dragen_anno=config["tools"]["cphi_dragen_anno"]
    shell:
        """
        (python3 {params.cphi_dragen_anno}/workflow/scripts/slivar/postfilter.py \
        --mode {wildcards.p} \
        --rare-main-vcf {input.rare_main} \
        --rare-clinvar-vcf {input.rare_clinvar} \
        --common-pathogenic-clinvar-vcf {input.common_pathogenic_clinvar} \
        --impact-order-file {params.cphi_dragen_anno}/workflow/scripts/slivar/default-order.txt \
        --out-vcf {output.vcf} &&
        bcftools sort -O v -o {output.vcf}.sorted {output.vcf} &&
        mv {output.vcf}.sorted {output.vcf}) > {log} 2>&1
        """


rule slivar_report:
    input:
        vcf="small_variants_slivar/{p}/{family}/{family}.{p}.postfilter.vcf"
    output:
        report=temp("small_variants_slivar/{p}/{family}/{family}.{p}.slivar.csv")
    wildcard_constraints:
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.report.log"
    conda:
        os.path.join(os.path.expanduser(config["tools"]["cphi_dragen_anno"]), "workflow", "envs", "slivar.yaml")
    params:
        cphi_dragen_anno=config["tools"]["cphi_dragen_anno"],
        crg2_pacbio=config["tools"]["crg2_pacbio"],
        hgmd=os.path.join(config["annotation"]["cre"]["database_path"], "hgmd_hg38.csv"),
    shell:
        """
        (python3 {params.cphi_dragen_anno}/workflow/scripts/slivar/build_report.py \
        --profile pacbio \
        --mode {wildcards.p} \
        --vcf {input.vcf} \
        --out-csv {output.report} \
        --impact-order-file {params.cphi_dragen_anno}/workflow/scripts/slivar/default-order.txt \
        --cre-data-dir {params.crg2_pacbio}/scripts/cre/data \
        --hgmd {params.hgmd}) > {log} 2>&1
        """


def get_cre_report_for_slivar_comparison(wildcards):
    report_type = {
        "coding": "wgs.coding",
        "wgs-high-impact": "wgs.high.impact",
        "panel": "panel",
        "panel-flank": "panel-flank",
    }[wildcards.p]
    suffix = ".SF" if acmg_sf_enabled and wildcards.p in {
        "coding", "wgs-high-impact"
    } else ""
    return f"reports/{wildcards.family}.{report_type}.CH{suffix}.csv"


def get_slivar_report_for_comparison(wildcards):
    report_type = {
        "coding": "wgs.coding",
        "wgs-high-impact": "wgs.high.impact",
        "panel": "panel",
        "panel-flank": "panel-flank",
    }[wildcards.p]
    suffix = ".SF" if acmg_sf_enabled and wildcards.p in {
        "coding", "wgs-high-impact"
    } else ""
    return f"reports_slivar/{wildcards.family}.{report_type}.CH{suffix}.csv"


rule compare_slivar_report_keys:
    input:
        gemini=get_cre_report_for_slivar_comparison,
        slivar=get_slivar_report_for_comparison,
    output:
        summary="reports_slivar_compare/{family}.{p}.summary.tsv",
        shared="reports_slivar_compare/{family}.{p}.shared.csv",
        gemini_only="reports_slivar_compare/{family}.{p}.gemini_only.csv",
        slivar_only="reports_slivar_compare/{family}.{p}.slivar_only.csv",
    wildcard_constraints:
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.compare.log"
    conda:
        os.path.join(os.path.expanduser(config["tools"]["cphi_dragen_anno"]), "workflow", "envs", "slivar.yaml")
    params:
        cphi_dragen_anno=config["tools"]["cphi_dragen_anno"],
        out_prefix="reports_slivar_compare/{family}.{p}",
    shell:
        """
        (python3 {params.cphi_dragen_anno}/workflow/scripts/slivar/compare_report_variant_keys.py \
        --gemini-report {input.gemini} \
        --slivar-report {input.slivar} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """
