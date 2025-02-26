include: "rules/common.smk"
include: "rules/methbat.smk"


samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)
case = samples[samples["case_or_control"] == "case"]
controls = samples[samples["case_or_control"] == "control"]

rule all:
    input:
       expand("methbat/outliers/{sample}.annotated.outliers.csv", sample=case.index)