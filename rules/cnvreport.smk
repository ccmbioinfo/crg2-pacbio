rule bcftools_merge:
    input: get_cnv_dir
    output:
        vcf = "cnv/{family}.cnv.vcf.gz"
    log:
        "logs/cnv/{family}.cnv.bcftools.merge.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        (
        n_vcfs=$(ls {input}/*hificnv*.vcf.gz | wc -l)
        if [ "$n_vcfs" -eq 1 ]; then
            cp {input}/*hificnv*.vcf.gz {output}
            tabix {output}
        else
            bcftools merge -m none {input}/*hificnv*.vcf.gz | bgzip > {output}
            tabix {output}
        fi
        ) > {log} 2>&1
        """

rule truvari_collapse:
    input: "cnv/{family}.cnv.vcf.gz"
    output:
        merged_variants = "cnv/{family}.cnv.truvari.merge.vcf",
        collapsed_variants = temp("cnv/{family}.cnv.truvari.collapse.vcf")
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/cnv/{family}.cnv.truvari.merge.log"
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        (truvari collapse -i {input} \
                -o {output.merged_variants} \
                -c {output.collapsed_variants} \
                -f {params.ref} \
                --sizemin 50 --sizemax 50000000 --refdist 2000 --pctsize 0.5 --pctovl 0.5 --pctseq 0 -k common)  > {log} 2>&1
        """

rule fix_hifi_cnv_CI:
    input: "cnv/{family}.cnv.truvari.merge.vcf"
    output: temp("cnv/{family}.cnv.truvari.merge.fix.CIPOS.vcf")
    log: "logs/cnv/{family}.fix.CIPOS.log"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
    conda:
        "../envs/str_sv.yaml"
    shell:
        """
        (python3 {params.crg2_pacbio}/scripts/fix_hifi_cnv_CI.py -input_vcf {input} -output_vcf {output}) > {log} 2>&1
        """

rule cnv_snpeff:
    input: "cnv/{family}.cnv.truvari.merge.fix.CIPOS.vcf"
    output:
        vcf = "cnv/{family}.cnv.snpeff.vcf",
    log:
        "logs/cnv/{family}.snpeff.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["data_dir"]
    wrapper:
        get_wrapper_path("snpeff")

rule cnv_annotsv:
    input: "cnv/{family}.cnv.truvari.merge.fix.CIPOS.vcf"
    output:
        annotsv_annotated =  "cnv/{family}.AnnotSV.tsv",
        annotsv_unannotated =  temp("cnv/{family}.AnnotSV.unannotated.tsv")
    log: "logs/cnv/{family}.annotsv.log"
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

rule cnv_report:
    input: 
        cnv_index = "cnv/{family}.cnv.truvari.merge.vcf.gz.tbi",
        snpeff = "cnv/{family}.cnv.snpeff.vcf",
        annotsv = "cnv/{family}.AnnotSV.tsv"
    output: "cnv/{family}.cnv.csv"
    log: "logs/cnv/{family}.cnv.report.log"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
        omim = config["annotation"]["omim_path"],
        exon = config["annotation"]["sv_report"]["exon"],
        anno_path = config["annotation"]["sv_report"]["anno_path"],
        repeats = config["trgt"]["adotto_repeats"],
        cnv_inhouse_c4r = config["annotation"]["sv_report"]["cnv_inhouse_c4r"],
        cnv_inhouse_tg = config["annotation"]["sv_report"]["cnv_inhouse_tg"],
        gnomad_SV = config["annotation"]["sv_report"]["gnomad_SV"],
        dgv = config["annotation"]["sv_report"]["dgv"],
        ensembl = config["trgt"]["ensembl"],
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
                        -variant_type CNV \
                        -omim {params.omim} \
                        -exon {params.exon} \
                        -gnomad {params.gnomad_SV} \
                        -dgv {params.dgv} \
                        -ensembl {params.ensembl} \
                        -cnv_inhouse_c4r {params.cnv_inhouse_c4r} \
                        -cnv_inhouse_tg {params.cnv_inhouse_tg} \
                        -colorsdb {params.colorsdb} \
                        -odd_regions {params.anno_path}/GRCh38.oddRegions.bed \
                        -repeats {params.anno_path}/human_GRCh38_no_alt_analysis_set.trgt.bed \
                        -c4r {params.c4r} \
                        -dark_regions {params.anno_path}/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed \
                        -clingen_HI {params.anno_path}/ClinGen_haploinsufficiency_gene_GRCh38.bed \
                        -clingen_TS {params.anno_path}/ClinGen_triplosensitivity_gene_GRCh38.bed \
                        -clingen_disease {params.anno_path}/ClinGen_tableExport_202310.csv \
                        -clingen_regions {params.anno_path}/ClinGen_region_curation_list_GRCh38.tsv) > {log} 2>&1
            else
                (python3 {params.crg2_pacbio}/scripts/annotate_SVs.py \
                    -annotsv {input.annotsv} \
                    -snpeff {input.snpeff} \
                    -variant_type CNV \
                    -omim {params.omim} \
                    -hpo {params.HPO} \
                    -exon {params.exon} \
                    -gnomad {params.gnomad_SV} \
                    -dgv {params.dgv} \
                    -ensembl {params.ensembl} \
                    -cnv_inhouse_c4r {params.cnv_inhouse_c4r} \
                    -cnv_inhouse_tg {params.cnv_inhouse_tg} \
                    -colorsdb {params.colorsdb} \
                    -odd_regions {params.anno_path}/GRCh38.oddRegions.bed \
                    -repeats {params.anno_path}/human_GRCh38_no_alt_analysis_set.trgt.bed \
                    -c4r {params.c4r} \
                    -dark_regions {params.anno_path}/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed \
                    -clingen_HI {params.anno_path}/ClinGen_haploinsufficiency_gene_GRCh38.bed \
                    -clingen_TS {params.anno_path}/ClinGen_triplosensitivity_gene_GRCh38.bed \
                    -clingen_disease {params.anno_path}/ClinGen_tableExport_202310.csv \
                    -clingen_regions {params.anno_path}/ClinGen_region_curation_list_GRCh38.tsv) > {log} 2>&1
            fi
        """
