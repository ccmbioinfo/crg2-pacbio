PIPELINE_VERSION="0.9.0"

include: "rules/common.smk"
include: "rules/annotation.smk"
include: "rules/snvreport.smk"

samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)

##### Target rules #####
project = config["run"]["project"]
family = project

rule all:
    input:
        "report/coding/{family}".format(family=project)