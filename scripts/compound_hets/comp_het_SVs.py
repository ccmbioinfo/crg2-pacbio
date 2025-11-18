import argparse
from functools import reduce
import pandas as pd
import numpy as np
import re

import compound_hets

def main():
    parser = argparse.ArgumentParser(
        description="Annotate compound het status for SVs."
    )
    parser.add_argument(
        "--sv", type=str, required=True, help="Path to SV file"
    )
    parser.add_argument(
        "--pedigree", type=str, required=True, help="Path to pedigree file"
    )
    args = parser.parse_args()

    # extract family and proband information
    family = args.pedigree.split("/")[-1].split(".")[0].split("_")[0]
    fam_dict = compound_hets.infer_pedigree_roles(args.pedigree)
    proband_id = fam_dict["child"]

    # read in SV file and filter
    SV = pd.read_csv(args.sv, low_memory=False)
    SV["Variant_id"] = SV["CHROM"] + "-" + SV["POS"].astype(str) + "-" + SV["END"].astype(str) +  "-" + SV["SVTYPE"] + "-" + SV["ID"]
    SV_rare_high_impact = SV[(SV["VARIANT"] != "intergenic_region") & (SV["gnomad_maxAF"] <= 0.01)].copy()
    SV_rare_high_impact.rename(columns={"ENSEMBL_GENE": "Ensembl_gene_id"}, inplace=True) # for compatibility with SNV CH functions
    # filter for variants that are het in the proband
    SV_rare_high_impact = SV_rare_high_impact[SV_rare_high_impact[f"{proband_id}_zyg"] == "het"]

    # get per-sample SV genotype status
    # extract all sample IDs based on genotype columns
    sample_ids = sorted(col.rsplit("_", 1)[0] for col in SV.columns if "_GT" in col)
    variant_ps = compound_hets.melt_sample_columns_SV(SV_rare_high_impact, "Variant_id", "_PS", "PS", sample_ids)
    variant_gts = compound_hets.melt_sample_columns_SV(
        SV_rare_high_impact, "Variant_id", "_GT", "GT", sample_ids
    )
    variant_zyg = compound_hets.melt_sample_columns_SV(
        SV_rare_high_impact, "Variant_id", "_zyg", "Zygosity", sample_ids
    )
    # combine all melted dataframes on Variant_id and Sample to create one dataframe of variants in long format
    dfs = [variant_ps, variant_gts, variant_zyg]
    # successively merges a list of DataFrames (dfs), each containing different variant/sample-level data,
    # to produce a single DataFrame (variant_gt_details) with all variant information for each (Variant_id, Sample) pair.
    variant_gt_details = reduce(
        lambda left, right: left.merge(
            right, on=["Variant_id", "Sample"], how="outer"
        ),
        dfs,
    )

    # create variant_to_gene table: variant ID and gene ID
    variant_to_gene = SV_rare_high_impact[["Variant_id", "Ensembl_gene_id"]].drop_duplicates()
    variant_to_gene = variant_to_gene[variant_to_gene['Ensembl_gene_id'].notna()]
    # NOTE:snpeff annotates some genes against transcript IDs or JASPAR matric IDs - should I filter these out?
    # an SV may overlap multiple genes, so we need to split the gene IDs and explode the dataframe to create one row per gene
    variant_to_gene["Ensembl_gene_id"] = variant_to_gene["Ensembl_gene_id"].str.split(";")
    variant_to_gene = variant_to_gene.explode("Ensembl_gene_id")
    variant_to_gene["Ensembl_gene_id"] = variant_to_gene["Ensembl_gene_id"].str.split("-") # some SVs annotated as gene-gene by SnpEff, e.g. if it's a feature fusion
    variant_to_gene = variant_to_gene.explode("Ensembl_gene_id")
    # now add gene IDs to variant_gt_details, which will have one row per variant per gene per sample
    variant_gt_details = variant_gt_details.merge(
        variant_to_gene, on="Variant_id", how="left"
    )
    # recode zygosity to be consistent with SNV CH functions
    variant_gt_details.replace({"het": "heterozygous", "hom": "homozygous_alt", "-": "homozygous_ref", "./.": "missing"}, inplace=True)
    variant_gt_details["GT_abstracted"] = variant_gt_details["GT"].replace({"0|1": "ref|alt", "1|0": "alt|ref"})
    # restrict variants to those that are het in the proband
    variant_gt_details_proband = variant_gt_details[
        (variant_gt_details["Sample"] == proband_id)
        & (variant_gt_details["Zygosity"] == "heterozygous")
    ]
    
    # determine gene compound het status using only long-read phasing information
    compound_het_status_no_parents = compound_hets.determine_compound_het_status_no_parents(
        variant_to_gene, variant_gt_details_proband
    )
    print(compound_het_status_no_parents["compound_het_status"].value_counts())
    # now if any gene has unknown compound het status, we attempt to determine the compound het status using parental genotypes
    # first, get genes with unknown compound het status
    unknown_gene = compound_het_status_no_parents[
        compound_het_status_no_parents["compound_het_status"] == "Unknown"
    ]
    unknown_variants = variant_gt_details[
        variant_gt_details["Ensembl_gene_id"].isin(unknown_gene["Ensembl_gene_id"])
    ].copy().drop_duplicates()

    # use parental genotypes to determine compound het status
    gene_haplotype_counts = compound_hets.count_gene_variants_by_parental_haplotype(
        unknown_variants, fam_dict
    )
    gene_haplotype_counts["compound_het_status"] = gene_haplotype_counts.apply(
        lambda x: compound_hets.is_compound_het(
            x["maternal_haplotype_count"],
            x["paternal_haplotype_count"],
            x["unknown_haplotype_count"],
        ),
        axis=1,
        )
    compound_het_status = compound_het_status_no_parents.copy().set_index(
        "Ensembl_gene_id"
    )
    for gene in gene_haplotype_counts.index:
        compound_het_status.loc[gene, "compound_het_status"] = (
            gene_haplotype_counts.loc[gene, "compound_het_status"]
        )
    compound_het_status.compound_het_status.value_counts()

    # export table of variants in genes with compound het status
    variant_gt_details = variant_gt_details.merge(compound_het_status, on="Ensembl_gene_id", how="left")
    variant_gt_details.to_csv(f"{family}_compound_het_SVs.csv", index=False)

main()

 