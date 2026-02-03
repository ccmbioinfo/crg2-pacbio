include: "rules/common.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"
include: "rules/svreport.smk"
include: "rules/outlier_expansions.smk"
include: "rules/pathogenic_expansion_loci.smk"
include: "rules/denovo_TR.smk"
include: "rules/compound_hets.smk"
include: "rules/cnvreport.smk"

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
if config["run"]["ped"]:
	children = get_children_ids(config["run"]["ped"])

rule all:
    input:
        "reports/{family}.wgs.coding.CH.csv".format(family=project),
        "reports/{family}.sv.CH.csv".format(family=project),
        "reports/{family}.cnv.CH.csv".format(family=project),
        "reports/{family}.compound.het.status.CH.csv".format(family=project),
        "reports/{family}.panel.CH.csv".format(family=project),
        "reports/{family}.panel-flank.CH.csv".format(family=project),
        "reports/{family}.wgs.high.impact.CH.csv".format(family=project),
        "reports/{family}.repeat.outliers.annotated.csv".format(family=project),
        "reports/{family}.known.path.str.loci.csv".format(family=project),
        expand("reports/{family}_{child}.TRGT.denovo.annotated.csv",
               family=project,
               child=children) if len(children) > 0 else []

