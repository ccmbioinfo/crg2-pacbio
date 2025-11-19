import argparse
from functools import reduce
import pandas as pd
import numpy as np
import re

import compound_hets

def main():
    parser = argparse.ArgumentParser(
        description="Annotate compound het status across variant types."
    )
    parser.add_argument(
        "--high_med", type=str, required=True, help="Path to HIGH-MED sequence variants TSV"
    )
    parser.add_argument(
        "--low", type=str, required=True, help="Path to LOW sequence variants TSV"
    )
    parser.add_argument(
        "--sequence_variant_report", type=str, required=True, help="Path to sequence variant report CSV"
    )
    parser.add_argument(
        "--sv", type=str, required=True, help="Path to SV report TSV"
    )
    parser.add_argument(
        "--pedigree", type=str, required=True, help="Path to pedigree file"
    )

    args = parser.parse_args()
    family = args.pedigree.split("/")[-1].split(".")[0].split("_")[0]
    fam_dict = compound_hets.infer_pedigree_roles(args.pedigree)
    proband_id = fam_dict["child"]

    # first process sequence variants
    high_med = pd.read_csv(args.high_med, sep="\t")
    low = pd.read_csv(args.low, sep="\t", low_memory=False)

    # process low impact variants to select genic variants with relatively high CADD or SpliceAI scores (or indels without these scores)
    low_impact_var_filter_scores = compound_hets.filter_low_impact_variants(low)
    # create high_impact table: genic damaging noncoding variants and nonsynonymous variants
    high_impact = pd.concat([low_impact_var_filter_scores, high_med])
    # de-duplicate (may be LOW impact variants in high_med that are also in low_impact if they have ClinVar annotations)
    high_impact = high_impact.drop_duplicates(subset=["Variant_id"])
    # create variant ID in the form chr-pos-ref-alt
    high_impact["Variant_id"] = high_impact["Chrom"] + "-" + high_impact["Pos"].astype(str) + "-" + high_impact["Ref"].astype(str) + "-" + high_impact["Alt"].astype(str)
    high_impact["Variant_id"] = high_impact["Variant_id"].apply(lambda x: x.replace("chr", ""))

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
    variant_qual = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gt_quals.", "GT_qual", sample_ids
    )

    # combine all melted dataframes on Variant_id and Sample to create one dataframe of variants in long format
    dfs = [variant_ps, variant_gts, variant_gt_type, variant_qual]
    # successively merges a list of DataFrames (dfs), each containing different variant/sample-level data,
    # to produce a single DataFrame (variant_gt_details) with all variant information for each (Variant_id, Sample) pair.
    sequence_variant_gt_details = reduce(
        lambda left, right: left.merge(
            right, on=["Variant_id", "Ref", "Alt", "Sample"], how="outer"
        ),
        dfs,
    )
    sequence_variant_gt_details.loc[sequence_variant_gt_details["GT"] == "./.", "GT_type"] = (
        2  # for some reason gemini codes some missing genotypes types as 3
    )
    sequence_variant_gt_details["Zygosity"] = sequence_variant_gt_details["GT_type"].map(
        lambda x: (
            "homozygous_alt"
            if x == 3
            else "heterozygous" if x == 1 else "homozygous_ref" if x == 0 else "missing"
        )
    )
    # create variant_to_gene table: variant ID and gene ID
    sequence_variant_to_gene = high_impact[["Variant_id", "Ensembl_gene_id"]].drop_duplicates()
    sequence_variant_gt_details = sequence_variant_gt_details.merge(
        sequence_variant_to_gene, on="Variant_id", how="left"
    )

    # abstract the genotype to be "ref|alt" or "alt|ref" to facilitate compound het status determination
    sequence_variant_gt_details["GT_abstracted"] = sequence_variant_gt_details.apply(
        lambda x: compound_hets.abstract_gt(x["Ref"], x["Alt"], x["GT"]), axis=1
    )

    # now process SVs
    # read in SV file and filter for rare genic variants
    SV = pd.read_csv(args.sv, low_memory=False)
    SV_rare_high_impact = SV[(SV["VARIANT"] != "intergenic_region") & (SV["gnomad_maxAF"] <= 0.01)].copy()
    SV_rare_high_impact["Variant_id"] = SV_rare_high_impact["CHROM"] + "-" + SV_rare_high_impact["POS"].astype(str) + "-" + SV_rare_high_impact["END"].astype(str) +  "-" + SV_rare_high_impact["SVTYPE"] + "-" + SV_rare_high_impact["ID"]
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
    SV_gt_details = reduce(
        lambda left, right: left.merge(
            right, on=["Variant_id", "Sample"], how="outer"
        ),
        dfs,
    )
    # create variant_to_gene table: variant ID and gene ID
    SV_to_gene = SV_rare_high_impact[["Variant_id", "Ensembl_gene_id"]].drop_duplicates()
    SV_to_gene = SV_to_gene[SV_to_gene['Ensembl_gene_id'].notna()]
    # NOTE:snpeff annotates some genes against transcript IDs or JASPAR matric IDs - should I filter these out?
    # an SV may overlap multiple genes, so we need to split the gene IDs and explode the dataframe to create one row per gene
    SV_to_gene["Ensembl_gene_id"] = SV_to_gene["Ensembl_gene_id"].str.split(";")
    SV_to_gene = SV_to_gene.explode("Ensembl_gene_id")
    SV_to_gene["Ensembl_gene_id"] = SV_to_gene["Ensembl_gene_id"].str.split("-") # some SVs annotated as gene-gene by SnpEff, e.g. if it's a feature fusion
    SV_to_gene = SV_to_gene.explode("Ensembl_gene_id")
    # now add gene IDs to variant_gt_details, which will have one row per variant per gene per sample
    SV_gt_details = SV_gt_details.merge(
        SV_to_gene, on="Variant_id", how="left"
    )
    # recode zygosity to be consistent with SNV CH functions
    SV_gt_details.replace({"het": "heterozygous", "hom": "homozygous_alt", "-": "homozygous_ref", "./.": "missing"}, inplace=True)
    SV_gt_details["GT_abstracted"] = SV_gt_details["GT"].replace({"0|1": "ref|alt", "1|0": "alt|ref"})


    # now combine SNV and SV dataframes and identify compound het status for each gene
    sequence_variant_gt_details = sequence_variant_gt_details[["Variant_id", "Sample", "PS", "GT_abstracted", "Zygosity", "Ensembl_gene_id"]]
    sequence_variant_gt_details["Variant_type"] = "sequence_variant"
    SV_gt_details = SV_gt_details[["Variant_id", "Sample", "PS", "GT_abstracted", "Zygosity", "Ensembl_gene_id"]]
    SV_gt_details["Variant_type"] = "SV"
    # combine all variant types into one dataframe
    all_variants = pd.concat([sequence_variant_gt_details, SV_gt_details])
    # remove variants without annotated gene 
    all_variants = all_variants[~all_variants["Ensembl_gene_id"].isna()]
    # get variants that are heterozygous in the proband
    all_variants_proband = all_variants[(all_variants["Sample"] == proband_id) & (all_variants["Zygosity"] == "heterozygous")].copy()
    # attempt to phase variants using proband phasing only
    print("Attempting to phase variants using proband phasing only for variants")
    compound_het_status_no_parents = compound_hets.determine_compound_het_status_no_parents(all_variants_proband)
    # now if any gene has unknown compound het status, we attempt to determine the compound het status using parental genotypes
    # first, get genes with unknown compound het status
    unknown_gene = compound_het_status_no_parents[
        compound_het_status_no_parents["compound_het_status"] == "UNKNOWN"
    ]
    unknown_variants = all_variants[
        all_variants["Ensembl_gene_id"].isin(unknown_gene["Ensembl_gene_id"])
    ].copy().drop_duplicates()
    # use parental genotypes to determine compound het status
    print("Attempting to phase variants in genes with unknown compound het status using parental genotypes")
    compound_het_status_parents = compound_hets.count_gene_variants_by_parental_haplotype(
        unknown_variants, fam_dict
    )
    compound_het_status_parents["compound_het_status"] = compound_het_status_parents.apply(
        lambda x: compound_hets.is_compound_het(
            x["maternal_haplotype_count"],
            x["paternal_haplotype_count"],
            x["unknown_haplotype_count"],
        ),
        axis=1,
        )
    # replace unknown compound het status with compound het status after parental phasing
    for gene in compound_het_status_parents.index:
        compound_het_status_no_parents.loc[compound_het_status_no_parents["Ensembl_gene_id"] == gene, "compound_het_status"] = (
            compound_het_status_parents.loc[gene, "compound_het_status"]
        )
    # add cross-variant compound het status to variants 
    all_variants = all_variants.merge(compound_het_status_no_parents, on="Ensembl_gene_id", how="left").sort_values(by=["Ensembl_gene_id", "Sample"])
    all_variants.to_csv(f"{family}_compound_het_variants.csv", index=False)
    
    gene_CH_status = all_variants[["Ensembl_gene_id", "Variant_type", "compound_het_status"]].drop_duplicates()
    gene_CH_status = gene_CH_status.groupby("Ensembl_gene_id").agg({"Variant_type": ";".join, "compound_het_status": "first"}).reset_index()
    gene_CH_status.rename(columns={"Variant_type": "CH_variant_types", "compound_het_status": "CH_status"}, inplace=True)

    # now annotate variant reports with compound het status 
    sequence_variant_report = pd.read_csv(args.sequence_variant_report)
    sequence_variant_report = sequence_variant_report.merge(gene_CH_status, on="Ensembl_gene_id", how="left")
    sequence_variant_report.fillna(".", inplace=True)
    sequence_variant_report.to_csv(f"{family}_sequence_variant_report_CH.csv", index=False)

    result_series = SV.apply(
        lambda x: compound_hets.SV_comp_het_status(x["ENSEMBL_GENE"], x["VARIANT"], gene_CH_status.set_index("Ensembl_gene_id")), 
        axis=1
    )
    SV[["CH_status", "CH_variant_types"]] = pd.DataFrame(result_series.tolist(), index=SV.index)
    SV.fillna(".", inplace=True)
    SV.to_csv(f"{family}_SV_report_CH.csv", index=False)



main()