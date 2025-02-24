include: "rules/common.smk"
include: "rules/methbat.smk"


samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)

rule all:
    input:
       expand("methbat/outliers/{sample}.annotated.outliers.csv", sample=samples.index)