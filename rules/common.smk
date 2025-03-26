import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version
from datetime import date

min_version("5.7.1")

report: "../report/workflow.rst"

###### Config file and sample sheets #####
#configfile: "config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)
case = samples[samples["case_or_control"] == "case"]
controls = samples[samples["case_or_control"] == "control"]
#validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["run"]["units"], dtype=str).set_index(["family"], drop=False)
#validate(units, schema="../schemas/units.schema.yaml")

variants_for_methbat = pd.read_table(config["run"]["variants_for_methbat"], dtype=str).set_index(["variant_type"], drop=False)

project = config["run"]["project"]

def get_wrapper_path(*dirs):
    return "file:%s" % os.path.join(workflow.basedir, "wrappers", *dirs)

def format_pedigree(wildcards):
    family = wildcards.family
    ped = config["run"]["ped"]
    if ped == "":
        return None
    else:
        ped = pd.read_csv(
            ped,
            sep=" ",
            header=None,
            names=["fam_id", "individual_id", "pat_id", "mat_id", "sex", "phenotype"],
        )

    for col in ["individual_id", "pat_id", "mat_id"]:
            ped[col] = [parse_ped_id(individual_id, family) for individual_id in ped[col].values]       
    
    ped["fam_id"] = family

    for col in ["individual_id", "pat_id", "mat_id"]:
        for row in range(len(ped[col])):
            if len(ped.loc[row, col].split("_")[-1]) != 1 and len(ped.loc[row, col].split("_")[-1]) != 2:
                pass   
            else:
                ped.loc[row, col] = ped.loc[row, col].replace(ped.loc[row, col], "0")

    ped = ped.drop(ped.index[row] for row in ped.index if ped["individual_id"][row] == '0' )


    ped.to_csv(f"{family}.ped", sep=" ", index=False, header=False)

    return f"{family}.ped"

def get_pbsv_vcf(wildcards):
    family = project
    input_vcf = units.loc[family, "pbsv_vcf"]

    return input_vcf

def get_trgt_path_str_vcf_dir(wildcards):
    family = project
    input_vcf = units.loc[family, "trgt_pathogenic_vcf_dir"]

    return input_vcf

def get_bam(wildcards):
    bam = samples.loc[wildcards.sample, "BAM"]

    return bam

def get_methbat_profiles(wildcards):
    return [f"methbat/profiles/{sample}.profile.tsv" for sample in controls.index]

# input functions for methbat outlier workflow
def get_sample_smvs(wildcards):
    return variants_for_methbat.loc["small_variants", "path"]

def get_pbsv_csv(wildcards):
    return variants_for_methbat.loc["SVs", "path"]

def get_cnvs(wildcards):
    return variants_for_methbat.loc["CNVs", "path"]

def get_TR_outliers (wildcards):
    return variants_for_methbat.loc["TR_outliers", "path"]