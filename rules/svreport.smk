rule snpeff:
    input: get_pbsv_vcf
    output:
        vcf = "sv/{family}.pbsv.snpeff.vcf",
    log:
        "logs/sv/{family}.snpeff.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["data_dir"]
    wrapper:
        get_wrapper_path("snpeff")

rule annotsv:
    input: get_pbsv_vcf
    output:
        annotsv_annotated =  "sv/{family}.AnnotSV.tsv",
        annotsv_unannotated =  temp("sv/{family}.AnnotSV.unannotated.tsv")
    log: "logs/sv/{family}.annotsv.log"
    params:
        annotsv_path = config["tools"]["annotSV"]
    conda:
        "../envs/common.yaml"
    shell:
        """
        (export ANNOTSV={params.annotsv_path}; 
        {params.annotsv_path}/bin/AnnotSV \
            -SVinputFile {input} \
            -outputFile {output.annotsv_annotated} \
            -overlap 50 \
            -reciprocal 1 \
            -genomeBuild GRCh38) > {log} 2>&1
        """

rule sv_report:
    input: 
        snpeff = "sv/{family}.pbsv.snpeff.vcf",
        annotsv = "sv/{family}.AnnotSV.tsv"
    output: "sv/{family}.sv.csv"
    log: "logs/sv/{family}.sv.report.log"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
        omim = config["annotation"]["omim_path"],
        exon = config["annotation"]["general"]["exon"],
        long_read_regions = config["annotation"]["general"]["long_read_regions"],
        repeats = config["annotation"]["general"]["adotto_repeats"],
        inhouse_c4r = config["annotation"]["sv_report"]["inhouse_c4r"],
        inhouse_tg = config["annotation"]["sv_report"]["inhouse_tg"],
        gnomad_SV = config["annotation"]["sv_report"]["gnomad_SV"],
        dgv = config["annotation"]["sv_report"]["dgv"],
        ensembl = config["annotation"]["general"]["ensembl"],
        clingen_path = config["annotation"]["general"]["clingen_path"],
        colorsdb = config["annotation"]["sv_report"]["colorsdb"],
        c4r = config["annotation"]["c4r"],
    conda:
        "../envs/str_sv.yaml"
    shell:
        """
        if [[ {params.HPO} == "none" ]]
            then
                    (python3 {params.crg2_pacbio}/scripts/annotate_SVs.py \
                        -annotsv {input.annotsv} \
                        -snpeff {input.snpeff} \
                        -variant_type SV \
                        -omim {params.omim} \
                        -exon {params.exon} \
                        -gnomad {params.gnomad_SV} \
                        -dgv {params.dgv} \
                        -ensembl {params.ensembl} \
                        -inhouse_c4r {params.inhouse_c4r} \
                        -inhouse_tg {params.inhouse_tg} \
                        -colorsdb {params.colorsdb} \
                        -odd_regions {params.long_read_regions}/GRCh38.oddRegions.bed \
                        -repeats {params.repeats} \
                        -c4r {params.c4r} \
                        -dark_regions {params.long_read_regions}/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed \
                        -clingen_HI {params.clingen_path}/ClinGen_haploinsufficiency_gene_GRCh38.bed \
                        -clingen_TS {params.clingen_path}/ClinGen_triplosensitivity_gene_GRCh38.bed \
                        -clingen_disease {params.clingen_path}/ClinGen_tableExport_202310.csv \
                        -clingen_regions {params.clingen_path}/ClinGen_region_curation_list_GRCh38.tsv) > {log} 2>&1
            else
                (python3 {params.crg2_pacbio}/scripts/annotate_SVs.py \
                    -annotsv {input.annotsv} \
                    -snpeff {input.snpeff} \
                    -variant_type SV \
                    -omim {params.omim} \
                    -hpo {params.HPO} \
                    -exon {params.exon} \
                    -gnomad {params.gnomad_SV} \
                    -dgv {params.dgv} \
                    -ensembl {params.ensembl} \
                    -inhouse_c4r {params.inhouse_c4r} \
                    -inhouse_tg {params.inhouse_tg} \
                    -colorsdb {params.colorsdb} \
                    -odd_regions {params.long_read_regions}/GRCh38.oddRegions.bed \
                    -repeats {params.repeats} \
                    -c4r {params.c4r} \
                    -dark_regions {params.long_read_regions}/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed \
                    -clingen_HI {params.clingen_path}/ClinGen_haploinsufficiency_gene_GRCh38.bed \
                    -clingen_TS {params.clingen_path}/ClinGen_triplosensitivity_gene_GRCh38.bed \
                    -clingen_disease {params.clingen_path}/ClinGen_tableExport_202310.csv \
                    -clingen_regions {params.clingen_path}/ClinGen_region_curation_list_GRCh38.tsv) > {log} 2>&1
            fi
        """
