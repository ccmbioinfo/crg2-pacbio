import argparse
from functools import reduce
import pandas as pd
import numpy as np
import re

import compound_hets

def main():
    parser = argparse.ArgumentParser(
        description="Annotate compound het status for sequence variants."
    )
    parser.add_argument(
        "--high_med", type=str, required=True, help="Path to HIGH-MED variants file"
    )
    parser.add_argument(
        "--low", type=str, required=True, help="Path to LOW variants file"
    )
    parser.add_argument(
        "--pedigree", type=str, required=True, help="Path to pedigree file"
    )

    args = parser.parse_args()

    family = args.pedigree.split("/")[-1].split(".")[0].split("_")[0]
    fam_dict = compound_hets.infer_pedigree_roles(args.pedigree)
    proband_id = fam_dict["child"]
    high_med = pd.read_csv(args.high_med, sep="\t")
    low = pd.read_csv(args.low, sep="\t", low_memory=False)

    # process low impact variants to select high impact variants
    low_impact_var_filter_scores = compound_hets.filter_low_impact_variants(low)
    # create high_impact table: genic damaging noncoding variants and nonsynonymous variants
    high_impact = pd.concat([low_impact_var_filter_scores, high_med])
    # de-duplicate (may be LOW impact variants in high_med that are also in low_impact if they have ClinVar annotations)
    high_impact["Variant_id"] = high_impact["Chrom"] + "-" + high_impact["Pos"].astype(str) + "-" + high_impact["Ref"].astype(str) + "-" + high_impact["Alt"].astype(str)
    high_impact["Variant_id"] = high_impact["Variant_id"].apply(lambda x: x.replace("chr", ""))
    high_impact = high_impact.drop_duplicates(subset=["Variant_id"])

    # restrict variants to those that are het in the proband
    high_impact = high_impact[high_impact[f"gt_types.{proband_id}"] == 1]
    # remove variants where genotype is like ./A etc
    high_impact = high_impact[~high_impact[f"gts.{proband_id}"].str.contains("\.")]

    # get per-sample variant genotype status
    # extract all sample IDs based on genotype columns
    sample_ids = sorted(
        {
            re.sub(r"^gts\.", "", col)
            for col in high_impact.columns
            if col.startswith("gts.")
        }
    )

    # split comma-separated PS string into individual PS columns for each sample
    ps_split_cols = high_impact["PS"].str.split(",", expand=True)
    for idx, sample in enumerate(sample_ids):
        high_impact[f"PS.{sample}"] = ps_split_cols[idx]
    high_impact.drop(columns=["PS"], inplace=True)

    variant_ps = compound_hets.melt_sample_columns(high_impact, "Variant_id", "PS.", "PS", sample_ids)
    variant_gts = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gts.", "GT", sample_ids
    )
    variant_gt_type = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gt_types.", "GT_type", sample_ids
    )
    variant_phase = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gt_phases.", "Phase", sample_ids
    )
    variant_qual = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gt_quals.", "GT_qual", sample_ids
    )

    # combine all melted dataframes on Variant_id and Sample to create one dataframe of variants in long format
    dfs = [variant_ps, variant_gts, variant_gt_type, variant_phase, variant_qual]
    # successively merges a list of DataFrames (dfs), each containing different variant/sample-level data,
    # to produce a single DataFrame (variant_gt_details) with all variant information for each (Variant_id, Sample) pair.
    variant_gt_details = reduce(
        lambda left, right: left.merge(
            right, on=["Variant_id", "Ref", "Alt", "Sample"], how="outer"
        ),
        dfs,
    )
    variant_gt_details.loc[variant_gt_details["GT"] == "./.", "GT_type"] = (
        2  # for some reason gemini codes some missing genotypes types as 3
    )
    variant_gt_details["Zygosity"] = variant_gt_details["GT_type"].map(
        lambda x: (
            "homozygous_alt"
            if x == 3
            else "heterozygous" if x == 1 else "homozygous_ref" if x == 0 else "missing"
        )
    )
    variant_gt_details.set_index(["Variant_id", "Sample"], inplace=True)
    # create variant_to_gene table: variant ID and gene ID
    variant_to_gene = high_impact[["Variant_id", "Ensembl_gene_id"]].drop_duplicates()
    variant_gt_details.reset_index(inplace=True)
    variant_gt_details = variant_gt_details.merge(
        variant_to_gene, on="Variant_id", how="left"
    )
 
    # abstract the genotype to be "ref|alt" or "alt|ref" to facilitate compound het status determination
    variant_gt_details["GT_abstracted"] = variant_gt_details.apply(
        lambda x: compound_hets.abstract_gt(x["Ref"], x["Alt"], x["GT"]), axis=1
    )

    # restrict variants to those that are het in the proband
    variant_gt_details_proband = variant_gt_details[
        variant_gt_details["Sample"] == proband_id
    ]
    variant_gt_details_proband = variant_gt_details_proband[
        variant_gt_details_proband["Zygosity"] == "heterozygous"
    ]

    # determine gene compound het status using only long-read phasing information
    compound_het_status_no_parents = compound_hets.determine_compound_het_status_no_parents(
        variant_to_gene, variant_gt_details_proband
    )
    print(
        f"CH gene status prior to parental phasing:\n {compound_het_status_no_parents.compound_het_status.value_counts()}"
    )
    # if any gene has unknown compound het status, we attempt to determine the compound het status using parental genotypes
    # first, get genes with unknown compound het status
    unknown_gene = compound_het_status_no_parents[
        compound_het_status_no_parents["compound_het_status"] == "Unknown"
    ]
    print(
        f"Number of genes with unknown compound het status after long-read phasing: {len(unknown_gene)}"
    )
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
    print(
        f"CH status of unknown genes after parental phasing:\n {gene_haplotype_counts['compound_het_status'].value_counts()}"
    )
    compound_het_status = compound_het_status_no_parents.copy().set_index(
        "Ensembl_gene_id"
    )
    for gene in gene_haplotype_counts.index:
        compound_het_status.loc[gene, "compound_het_status"] = (
            gene_haplotype_counts.loc[gene, "compound_het_status"]
        )
    compound_het_status.compound_het_status.value_counts()
    print(
        f"CH gene status after parental phasing:\n {compound_het_status.compound_het_status.value_counts()}"
    )

    # export table of variants in genes with compound het status
    # ch_variants_df = get_compound_het_variants(
    #     high_impact, variant_to_gene, compound_het_status, proband_id, sample_ids
    # )
    # ch_variants_df.to_csv(f"{family}_compound_het_variants.csv", index=False)
    variant_gt_details = variant_gt_details.merge(compound_het_status, on="Ensembl_gene_id", how="left")
    variant_gt_details.drop(columns=["GT_type", "Phase", "GT_qual"], inplace=True)
    variant_gt_details.to_csv(f"{family}_compound_het_sequence_variants.csv", index=False)

main()
