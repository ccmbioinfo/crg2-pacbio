PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"
include: "rules/svreport.smk"
include: "rules/outlier_expansions.smk"
include: "rules/pathogenic_expansion_loci.smk"
include: "rules/denovo_TR.smk"

samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)

##### Target rules #####
project = config["run"]["project"]
family = project

rule all:
    input:
        "small_variants/coding/{family}".format(family=project),
        "small_variants/panel/{family}".format(family=project),
        "small_variants/panel-flank/{family}".format(family=project),
        "sv/{family}.pbsv.csv".format(family=project),
        "repeat_outliers/{family}.repeat.outliers.annotated.csv".format(family=project),
        "pathogenic_repeats/{family}.known.path.str.loci.csv".format(family=project),
        "TRGT_denovo/{family}.TRGT.denovo.annotated.csv".format(family=project)

