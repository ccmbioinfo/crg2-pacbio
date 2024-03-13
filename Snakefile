rule all: 
    input:
        "repeat_outliers/annotated_repeat_outliers.csv"

include: "rules/outlier_expansions.smk"