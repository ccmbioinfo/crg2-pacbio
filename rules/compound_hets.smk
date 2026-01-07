rule get_sequence_variants_for_CH:
    input:
        gemini_db="annotated/coding/{family}-gemini.db"
    output:
        variants="small_variants/{family}.{severity}.impact.variants.tsv",
    params:
        severity="{severity}",
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    log:
        "logs/compound_hets/{family}.get.sequence.variants.for.CH.{severity}.log",
    conda:
        "../envs/gemini.yaml"
    wildcard_constraints:
        severity="HIGH-MED|LOW"
    shell:
        "{params.crg2_pacbio}/scripts/compound_hets/get_sequence_var_for_CH.sh {input.gemini_db} {params.severity} > {output.variants}"

rule identify_compound_hets:
    input:
        high_med_variants="small_variants/{family}.HIGH-MED.impact.variants.tsv",
        low_variants="small_variants/{family}.LOW.impact.variants.tsv",
        small_variant_report_dir="small_variants/coding/{family}",
        SV_report="sv/{family}.sv.csv",
        CNV_report="cnv/{family}.cnv.csv",
        ensembl=config["annotation"]["general"]["ensembl"],
        ensembl_to_NCBI_df=config["annotation"]["ensembl_to_NCBI_df"],
        HPO=config["run"]["hpo"],
        pedigree=config["run"]["ped"],
    output:
        sequence_variant_report_CH="reports/{family}.wgs.coding.CH.csv",
        SV_report_CH="reports/{family}.sv.CH.csv",
        CNV_report_CH="reports/{family}.cnv.CH.csv",
        compound_het_status="reports/{family}.compound.het.status.csv",
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"]
    conda:
        "../envs/str_sv.yaml"
    log:
        "logs/compound_hets/{family}.identify.compound.hets.log",
    shell:
        """
        (python3 {params.crg2_pacbio}/scripts/annotate_compound_hets.py --high_med {input.high_med_variants} \
        --low {input.low_variants} \
        --sv {input.SV_report}  \
        --cnv {input.CNV_report}  \
        --ensembl {input.ensembl}  \
        --ensembl_to_NCBI_df {input.ensembl_to_NCBI_df}  \
        --pedigree {input.pedigree}  \
        --sequence_variant_report_dir {input.small_variant_report_dir}  \
        --family {wildcards.family}) > {log} 2>&1
        """
 