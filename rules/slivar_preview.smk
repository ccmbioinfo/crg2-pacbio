DRAGEN_REPO = os.path.abspath(os.path.expanduser(config["tools"]["cphi_dragen_anno"]))
SHARED_SLIVAR_DIR = os.path.join(DRAGEN_REPO, "workflow", "scripts", "slivar")
SHARED_SLIVAR_ENV = os.path.join(DRAGEN_REPO, "workflow", "envs", "slivar.yaml")
SHARED_SLIVAR_WRAPPER = "file:" + os.path.join(
    DRAGEN_REPO, "workflow", "wrappers", "slivar"
)
PACBIO_ROOT = os.path.abspath(workflow.basedir)


rule slivar_select:
    input:
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz"
    output:
        rare_main="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz",
        rare_main_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz.tbi",
        rare_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz",
        rare_clinvar_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz.tbi",
        common_pathogenic_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz",
        common_pathogenic_clinvar_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz.tbi",
    wildcard_constraints:
        family=project,
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.select.log"
    params:
        js=os.path.join(SHARED_SLIVAR_DIR, "slivar_functions.js"),
        consequence_order_file=os.path.join(SHARED_SLIVAR_DIR, "default-order.txt"),
        mode="{p}",
        profile="pacbio",
    wrapper:
        SHARED_SLIVAR_WRAPPER


rule slivar_postfilter:
    input:
        rare_main="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz",
        rare_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz",
        common_pathogenic_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz",
    output:
        vcf="small_variants_slivar/{p}/{family}/{family}.{p}.postfilter.vcf",
    wildcard_constraints:
        family=project,
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.postfilter.log"
    conda:
        SHARED_SLIVAR_ENV
    shell:
        """
        (python3 {SHARED_SLIVAR_DIR}/postfilter.py \
        --mode {wildcards.p} \
        --rare-main-vcf {input.rare_main} \
        --rare-clinvar-vcf {input.rare_clinvar} \
        --common-pathogenic-clinvar-vcf {input.common_pathogenic_clinvar} \
        --impact-order-file {SHARED_SLIVAR_DIR}/default-order.txt \
        --out-vcf {output.vcf} &&
        bcftools sort -O v -o {output.vcf}.sorted {output.vcf} &&
        mv {output.vcf}.sorted {output.vcf}) > {log} 2>&1
        """


rule slivar_report:
    input:
        vcf="small_variants_slivar/{p}/{family}/{family}.{p}.postfilter.vcf"
    output:
        report="reports_slivar_raw/{family}.{p}.slivar.csv"
    wildcard_constraints:
        family=project,
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.report.log"
    conda:
        SHARED_SLIVAR_ENV
    params:
        hgmd=os.path.join(config["annotation"]["cre"]["database_path"], "hgmd_hg38.csv"),
        cre_data_dir=os.path.join(PACBIO_ROOT, "scripts", "cre", "data"),
    shell:
        """
        (mkdir -p reports_slivar_raw &&
        python3 {SHARED_SLIVAR_DIR}/build_report.py \
        --profile pacbio \
        --mode {wildcards.p} \
        --vcf {input.vcf} \
        --out-csv {output.report} \
        --impact-order-file {SHARED_SLIVAR_DIR}/default-order.txt \
        --cre-data-dir {params.cre_data_dir} \
        --hgmd {params.hgmd}) > {log} 2>&1
        """


rule slivar_select_compound_het_candidates:
    input:
        vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz"
    output:
        candidates="small_variants_slivar/compound-hets/{family}/{family}.compound-het.candidates.vcf.gz",
    wildcard_constraints:
        family=project,
    log:
        "logs/compound_hets/{family}.slivar.select.compound.het.candidates.log"
    params:
        js=os.path.join(SHARED_SLIVAR_DIR, "slivar_functions.js"),
        consequence_order_file=os.path.join(SHARED_SLIVAR_DIR, "default-order.txt"),
        mode="compound-hets",
        profile="pacbio",
    wrapper:
        SHARED_SLIVAR_WRAPPER


rule get_slivar_small_variants_for_CH:
    input:
        vcf="small_variants_slivar/compound-hets/{family}/{family}.compound-het.candidates.vcf.gz",
    output:
        high_med="small_variants_slivar/{family}.HIGH-MED.impact.variants.tsv",
        low="small_variants_slivar/{family}.LOW.impact.variants.tsv",
    wildcard_constraints:
        family=project,
    log:
        "logs/compound_hets/{family}.slivar.get.sequence.variants.for.CH.log"
    conda:
        SHARED_SLIVAR_ENV
    shell:
        """
        (mkdir -p small_variants_slivar &&
        python3 {SHARED_SLIVAR_DIR}/build_ch_tsv.py \
        --profile pacbio \
        --vcf {input.vcf} \
        --impact-order-file {SHARED_SLIVAR_DIR}/default-order.txt \
        --high-med-out {output.high_med} \
        --low-out {output.low}) > {log} 2>&1
        """


slivar_hpo_available = config["run"].get("hpo", "")

slivar_hpo_panel_inputs = {
    "panel_variant_report": "reports_slivar_raw/{family}.panel.slivar.csv",
    "panel_flank_variant_report": "reports_slivar_raw/{family}.panel-flank.slivar.csv",
    "HPO": config["run"]["hpo"],
} if slivar_hpo_available else {}

slivar_hpo_panel_outputs = {
    "panel_variant_report_CH": "reports_slivar/{family}.panel.CH.csv",
    "panel_flank_variant_report_CH": "reports_slivar/{family}.panel-flank.CH.csv",
} if slivar_hpo_available else {}


def get_slivar_hpo_panel_args(wildcards, input):
    if slivar_hpo_available:
        return (
            '--hpo "$hpo" '
            "--panel_variant_report_dir input/panel "
            "--panel_flank_variant_report_dir input/panel-flank"
        )
    return ""


def stage_slivar_hpo_panel_reports(wildcards, input):
    if slivar_hpo_available:
        return (
            f'hpo=$(realpath {input.HPO})\n'
            'mkdir -p "$stage/input/panel" "$stage/input/panel-flank"\n'
            f'ln -sf "$(realpath {input.panel_variant_report})" "$stage/input/panel/{wildcards.family}.panel.wgs.slivar.csv"\n'
            f'ln -sf "$(realpath {input.panel_flank_variant_report})" "$stage/input/panel-flank/{wildcards.family}.panel-flank.wgs.slivar.csv"'
        )
    return ""


def copy_slivar_hpo_panel_reports(wildcards):
    if slivar_hpo_available:
        return (
            f'cp -L reports/{wildcards.family}.panel.CH.csv "$root/reports_slivar/{wildcards.family}.panel.CH.csv"\n'
            f'cp -L reports/{wildcards.family}.panel-flank.CH.csv "$root/reports_slivar/{wildcards.family}.panel-flank.CH.csv"'
        )
    return ""


rule identify_compound_hets_slivar:
    input:
        high_med_variants="small_variants_slivar/{family}.HIGH-MED.impact.variants.tsv",
        low_variants="small_variants_slivar/{family}.LOW.impact.variants.tsv",
        small_variant_report="reports_slivar_raw/{family}.coding.slivar.csv",
        wgs_high_impact_variant_report="reports_slivar_raw/{family}.wgs-high-impact.slivar.csv",
        SV_report="sv/{family}.sv.csv",
        CNV_report="cnv/{family}.cnv.csv",
        ensembl=config["annotation"]["general"]["ensembl"],
        ensembl_to_NCBI_df=config["annotation"]["ensembl_to_NCBI_df"],
        pedigree=config["run"]["ped"],
        sample_order="small_variants/{family}.sample.order.txt",
        **slivar_hpo_panel_inputs,
    output:
        small_variant_report_CH=output_status("reports_slivar/{family}.wgs.coding.CH.csv"),
        wgs_high_impact_variant_report_CH=output_status("reports_slivar/{family}.wgs.high.impact.CH.csv"),
        SV_report_CH=output_status("reports_slivar/{family}.sv.CH.csv"),
        CNV_report_CH=output_status("reports_slivar/{family}.cnv.CH.csv"),
        compound_het_status="reports_slivar/{family}.compound.het.status.CH.csv",
        **slivar_hpo_panel_outputs,
    wildcard_constraints:
        family=project,
    params:
        seq_type="long",
        hpo_panel_args=get_slivar_hpo_panel_args,
        stage_hpo_panel_reports=stage_slivar_hpo_panel_reports,
        copy_hpo_panel_reports=copy_slivar_hpo_panel_reports,
        acmg_sf_flag=str(config["run"].get("acmg_sf", "false")).lower(),
        mavedb_tsv=config["annotation"]["general"]["mavedb_tsv"],
        annotate_compound_hets=os.path.join(PACBIO_ROOT, "scripts", "annotate_compound_hets.py"),
        add_mavedb=os.path.join(PACBIO_ROOT, "scripts", "add_mavedb_columns.py"),
    conda:
        "../envs/str_sv.yaml"
    log:
        "logs/compound_hets/{family}.slivar.identify.compound.hets.log"
    shell:
        """
        (set -e
        root=$(pwd)
        stage="$root/slivar_ch_work/{wildcards.family}"
        high_med=$(realpath {input.high_med_variants})
        low=$(realpath {input.low_variants})
        sample_order=$(realpath {input.sample_order})
        ensembl=$(realpath {input.ensembl})
        ensembl_to_ncbi=$(realpath {input.ensembl_to_NCBI_df})
        pedigree=$(realpath {input.pedigree})
        mavedb_tsv=$(realpath {params.mavedb_tsv})
        small_report=$(realpath {input.small_variant_report})
        high_impact_report=$(realpath {input.wgs_high_impact_variant_report})
        sv_report=$(realpath {input.SV_report})
        cnv_report=$(realpath {input.CNV_report})
        mkdir -p "$stage/input/wgs-coding" "$stage/input/wgs-high-impact" "$stage/reports" "$root/reports_slivar"
        ln -sf "$small_report" "$stage/input/wgs-coding/{wildcards.family}.wgs.coding.slivar.csv"
        ln -sf "$high_impact_report" "$stage/input/wgs-high-impact/{wildcards.family}.wgs.high.impact.slivar.csv"
        {params.stage_hpo_panel_reports}
        cd "$stage"
        python3 {params.annotate_compound_hets} \
        --seq_type {params.seq_type} \
        --high_med "$high_med" \
        --low "$low" \
        --sv "$sv_report" \
        --cnv "$cnv_report" \
        --ensembl "$ensembl" \
        --ensembl_to_NCBI_df "$ensembl_to_ncbi" \
        --pedigree "$pedigree" \
        {params.hpo_panel_args} \
        --sequence_variant_report_dir input/wgs-coding \
        --wgs_high_impact_variant_report_dir input/wgs-high-impact \
        --sample_order "$sample_order" \
        --family {wildcards.family} \
        --acmg_sf {params.acmg_sf_flag}
        cp -L reports/{wildcards.family}.wgs.coding.CH.csv "$root/{output.small_variant_report_CH}"
        cp -L reports/{wildcards.family}.wgs.high.impact.CH.csv "$root/{output.wgs_high_impact_variant_report_CH}"
        cp -L reports/{wildcards.family}.sv.CH.csv "$root/{output.SV_report_CH}"
        cp -L reports/{wildcards.family}.cnv.CH.csv "$root/{output.CNV_report_CH}"
        cp -L reports/{wildcards.family}.compound.het.status.CH.csv "$root/{output.compound_het_status}"
        {params.copy_hpo_panel_reports}
        python3 {params.add_mavedb} \
        --family {wildcards.family} \
        --reports-dir "$root/reports_slivar" \
        --suffix csv \
        --mavedb-tsv "$mavedb_tsv") > {log} 2>&1
        """


def get_legacy_report_for_slivar_comparison(wildcards):
    report_type = {
        "coding": "wgs.coding",
        "wgs-high-impact": "wgs.high.impact",
        "panel": "panel",
        "panel-flank": "panel-flank",
    }[wildcards.p]
    suffix = sf_suffix if wildcards.p in {"coding", "wgs-high-impact"} else ""
    return f"reports/{wildcards.family}.{report_type}.CH{suffix}.csv"


def get_slivar_report_for_comparison(wildcards):
    report_type = {
        "coding": "wgs.coding",
        "wgs-high-impact": "wgs.high.impact",
        "panel": "panel",
        "panel-flank": "panel-flank",
    }[wildcards.p]
    suffix = sf_suffix if wildcards.p in {"coding", "wgs-high-impact"} else ""
    return f"reports_slivar/{wildcards.family}.{report_type}.CH{suffix}.csv"


rule compare_slivar_report_keys:
    input:
        gemini=get_legacy_report_for_slivar_comparison,
        slivar=get_slivar_report_for_comparison,
    output:
        summary="reports_slivar_compare/{family}.{p}.summary.tsv",
        shared="reports_slivar_compare/{family}.{p}.shared.csv",
        gemini_only="reports_slivar_compare/{family}.{p}.gemini_only.csv",
        slivar_only="reports_slivar_compare/{family}.{p}.slivar_only.csv",
    wildcard_constraints:
        family=project,
        p="coding|wgs-high-impact|panel|panel-flank",
    log:
        "logs/slivar/{family}.{p}.compare.log"
    conda:
        SHARED_SLIVAR_ENV
    params:
        out_prefix="reports_slivar_compare/{family}.{p}",
    shell:
        """
        (python3 {SHARED_SLIVAR_DIR}/compare_report_variant_keys.py \
        --gemini-report {input.gemini} \
        --slivar-report {input.slivar} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


slivar_acmg_input_report_type = [
    "wgs.coding.CH",
    "wgs.high.impact.CH",
    "sv.CH",
    "cnv.CH",
]


rule add_acmg_sf_columns_slivar:
    input:
        report="reports_slivar/{family}.{input_report_type}.csv",
        acmg_sf_list=config["annotation"]["general"]["acmg_sf_list"],
    output:
        report="reports_slivar/{family}.{input_report_type}.SF.csv",
    wildcard_constraints:
        family=project,
        input_report_type="|".join([t.replace(".", r"\.") for t in slivar_acmg_input_report_type]),
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
        seq_type="long",
    log:
        "logs/report/acmg_sf/slivar/{family}.{input_report_type}.SF.log"
    conda:
        "../envs/str_sv.yaml"
    script:
        os.path.join(DRAGEN_REPO, "workflow", "scripts", "add_acmg_sf_columns.py")


rule create_acmg_sf_report_slivar:
    input:
        reports=lambda wildcards: expand(
            "reports_slivar/{family}.{input_report_type}.SF.csv",
            family=wildcards.family,
            input_report_type=slivar_acmg_input_report_type,
        ),
    output:
        report="reports_slivar/{family}.ACMG.SF.csv",
    wildcard_constraints:
        family=project,
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
    log:
        "logs/report/acmg_sf/slivar/{family}.acmg_sf_report.log"
    conda:
        "../envs/str_sv.yaml"
    script:
        "../scripts/create_acmg_sf_report.py"
