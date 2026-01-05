PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"
include: "rules/svreport.smk"
include: "rules/outlier_expansions.smk"
include: "rules/pathogenic_expansion_loci.smk"
include: "rules/denovo_TR.smk"
include: "rules/compound_hets.smk"

def get_children_ids(ped_file):
    import pandas as pd
    pedigree = pd.read_csv(ped_file, sep=" ", header=None, 
                          names=["family_ID", "individual_ID", "paternal_ID", "maternal_ID", "sex", "phenotype"])
    children = pedigree[pedigree["paternal_ID"] != "0"][pedigree["maternal_ID"] != "0"]["individual_ID"].apply(lambda x: x.split("_")[1]).values
    family = pedigree["family_ID"].iloc[0]

    return children

samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)

##### Target rules #####
project = config["run"]["project"]
if config["run"]["ped"]:
	children = get_children_ids(config["run"]["ped"])

rule all:
    input:
        "small_variants/coding/{family}".format(family=project),
        "small_variants/panel/{family}".format(family=project),
        "small_variants/panel-flank/{family}".format(family=project),
        "small_variants/wgs-high-impact/{family}".format(family=project),
        "sv/{family}.pbsv.csv".format(family=project),
        "repeat_outliers/{family}.repeat.outliers.annotated.csv".format(family=project),
        "pathogenic_repeats/{family}.known.path.str.loci.csv".format(family=project),
        expand("TRGT_denovo/{family}_{child}.TRGT.denovo.annotated.csv",
               family=project,
               child=children) if config["run"]["ped"] else []

