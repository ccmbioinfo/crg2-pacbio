include: "rules/common.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"
include: "rules/svreport.smk"
include: "rules/annotate_mt_var.smk"
include: "rules/outlier_expansions.smk"
include: "rules/pathogenic_expansion_loci.smk"
include: "rules/denovo_TR.smk"
include: "rules/compound_hets.smk"
include: "rules/cnvreport.smk"
include: "rules/qc.smk"

def get_children_ids(ped_file):
    import pandas as pd
    pedigree = pd.read_csv(ped_file, sep=" ", header=None, 
                          names=["family_ID", "individual_ID", "paternal_ID", "maternal_ID", "sex", "phenotype"])
    pedigree = pedigree.astype(str)
    children = pedigree[pedigree["paternal_ID"] != "0"][pedigree["maternal_ID"] != "0"]["individual_ID"].apply(lambda x: x.split("_")[1]).values
    family = pedigree["family_ID"].iloc[0]

    return children

samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)

    
##### Target rules #####
project = config["run"]["project"]
children = []
if config["run"]["ped"]:
    children = get_children_ids(config["run"]["ped"])

#conditionally include acmf_sf rule
if str(config["run"].get("acmg_sf", "")).lower() == "true":
    include: "rules/acmg_sf.smk"
    
acmg_sf_enabled = str(config["run"].get("acmg_sf", "")).lower() == "true"
sf_suffix = ".SF" if acmg_sf_enabled else ""
slivar_preview_enabled = str(config["run"].get("slivar_preview", "")).lower() == "true"

acmg_sf_report_output = ["reports/{family}.ACMG.SF.csv".format(family=project)] if acmg_sf_enabled else []

slivar_preview_outputs = [
    "reports_slivar/{family}.wgs.coding.CH{sf}.csv".format(family=project, sf=sf_suffix),
    "reports_slivar/{family}.wgs.high.impact.CH{sf}.csv".format(family=project, sf=sf_suffix),
    "reports_slivar/{family}.sv.CH{sf}.csv".format(family=project, sf=sf_suffix),
    "reports_slivar/{family}.cnv.CH{sf}.csv".format(family=project, sf=sf_suffix),
    "reports_slivar/{family}.compound.het.status.CH.csv".format(family=project),
    "reports_slivar_compare/{family}.coding.summary.tsv".format(family=project),
    "reports_slivar_compare/{family}.wgs-high-impact.summary.tsv".format(family=project),
    "reports_slivar/{family}.repeat.outliers.annotated.csv".format(family=project),
    "reports_slivar/{family}.known.path.str.loci.csv".format(family=project),
    "reports_slivar/{family}.multiqc_report.html".format(family=project),
    "reports_slivar/{family}.mito.csv".format(family=project),
] if slivar_preview_enabled else []

if slivar_preview_enabled and acmg_sf_enabled:
    slivar_preview_outputs.append(
        "reports_slivar/{family}.ACMG.SF.csv".format(family=project)
    )

if slivar_preview_enabled and len(children) > 0:
    slivar_preview_outputs.extend(
        expand("reports_slivar/{family}_{child}.TRGT.denovo.annotated.csv",
               family=project,
               child=children)
    )

hpo_reports = []
if config["run"].get("hpo", ""):
    hpo_reports = [
        "reports/{family}.panel.CH.csv".format(family=project),
        "reports/{family}.panel-flank.CH.csv".format(family=project),
    ]
    if slivar_preview_enabled:
        hpo_reports.extend([
            "reports_slivar/{family}.panel.CH.csv".format(family=project),
            "reports_slivar/{family}.panel-flank.CH.csv".format(family=project),
            "reports_slivar_compare/{family}.panel.summary.tsv".format(family=project),
            "reports_slivar_compare/{family}.panel-flank.summary.tsv".format(family=project),
        ])
    
rule all:
    input:
        "reports/{family}.wgs.coding.CH{sf}.csv".format(family=project, sf=sf_suffix),
        "reports/{family}.sv.CH{sf}.csv".format(family=project, sf=sf_suffix),
        "reports/{family}.cnv.CH{sf}.csv".format(family=project, sf=sf_suffix),
        "reports/{family}.compound.het.status.CH.csv".format(family=project, sf=sf_suffix),
        "reports/{family}.wgs.high.impact.CH{sf}.csv".format(family=project, sf=sf_suffix),
        "reports/{family}.repeat.outliers.annotated.csv".format(family=project),
        "reports/{family}.known.path.str.loci.csv".format(family=project),
        "reports/{family}.multiqc_report.html".format(family=project),
        *hpo_reports,
        "reports/{family}.mito.csv".format(family=project),
        *acmg_sf_report_output,
        *slivar_preview_outputs,
        expand("reports/{family}_{child}.TRGT.denovo.annotated.csv",
               family=project,
               child=children) if len(children) > 0 else []
