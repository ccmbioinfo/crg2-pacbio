import argparse
from functools import reduce
import pandas as pd
import numpy as np
import re

def infer_pedigree_roles(pedigree: str) -> dict:
    """Map pedigree roles to individual IDs."""
    pedigree = pd.read_csv(
        pedigree,
        sep=" ",
        header=None,
        names=[
            "family_ID",
            "individual_ID",
            "paternal_ID",
            "maternal_ID",
            "sex",
            "phenotype",
        ],
    )
    # just take ID of first child if there are multiple children
    child = child = pedigree[
        (pedigree["paternal_ID"] != "0") & (pedigree["maternal_ID"] != "0")
    ]["individual_ID"].values[0]
    father = pedigree[
        pedigree["individual_ID"]
        == pedigree[pedigree["individual_ID"] == child]["paternal_ID"].values[0]
    ]["individual_ID"].values[0]
    mother = pedigree[
        pedigree["individual_ID"]
        == pedigree[pedigree["individual_ID"] == child]["maternal_ID"].values[0]
    ]["individual_ID"].values[0]
    fam_dict = {"child": child, "father": father, "mother": mother}

    return fam_dict

def determine_compound_het_status_no_parents(
    variant_to_gene: pd.DataFrame, variant_gt_details: pd.DataFrame
) -> pd.DataFrame:
    """Determine compound het status for a given gene using only long-read phasing information.

    Args:
        variant_to_gene (pd.DataFrame): DataFrame mapping variant IDs to Ensembl gene IDs
        variant_gt_details (pd.DataFrame): DataFrame with variant genotype details

    Returns:
        pd.DataFrame: DataFrame with compound het status for each gene
    """
    gene_to_compound_hets = {}
    for gene in variant_to_gene["Ensembl_gene_id"].unique():
        variants = variant_gt_details[variant_gt_details["Ensembl_gene_id"] == gene]
        num_variants = len(variants)  # total number of het variants in the gene
        true_compound_hets = 0
        unknown_compound_hets = 0
        phase_sets = variants["PS"].unique()
        if num_variants > 1:
            if phase_sets.tolist() == ["."]:  # variants are not phased
                unknown_compound_hets += 1
            elif len(phase_sets) == 1:  # all variants are in the same phase set
                PS = phase_sets[0]
                variants_in_PS_gts = variants[variants["PS"] == PS][
                    "GT_abstracted"
                ].unique()
                if "ref|alt" in variants_in_PS_gts and "alt|ref" in variants_in_PS_gts:
                    true_compound_hets += 1
            else:
                # if two phase sets, and one is unphased, can be compound het (if phased variant(s) are in the same phase set) or unknown
                # if two or more phase sets, can be compound het (if phased variant(s) are in the same phase set) or unknown
                unknown_compound_hets += 1  # if there's an unphased variant, we don't know if it's a compound het with phased variant(s)
                for PS in variants["PS"].unique():
                    if PS != ".":
                        variants_in_PS_gts = variants[variants["PS"] == PS][
                            "GT_abstracted"
                        ].unique()
                        if (
                            "ref|alt" in variants_in_PS_gts
                            and "alt|ref" in variants_in_PS_gts
                        ):
                            true_compound_hets += 1

        if true_compound_hets > 0:
            compound_het_status = "True"
        else:
            if unknown_compound_hets > 0:
                compound_het_status = "Unknown"
            else:
                compound_het_status = "False"
        gene_to_compound_hets[gene] = compound_het_status
    gene_to_compound_hets = pd.DataFrame.from_dict(
        gene_to_compound_hets, orient="index"
    ).reset_index()
    gene_to_compound_hets.columns = ["Ensembl_gene_id", "compound_het_status"]

    return gene_to_compound_hets


def melt_sample_columns(
    df: pd.DataFrame, id_col: str, prefix: str, value_name: str, sample_ids: list
) -> pd.DataFrame:
    """Generate melted/long format dataframefor a given sample-specific field's prefix (e.g. genotype), for all samples detected"""
    sample_cols = [f"{prefix}{sample}" for sample in sample_ids]
    melted = df.melt(
        id_vars=[id_col, "Ref", "Alt"],
        value_vars=sample_cols,
        var_name="Sample",
        value_name=value_name,
    )
    melted["Sample"] = melted["Sample"].str.replace(f"{prefix}", "", regex=False)

    return melted


def melt_sample_columns_SV(
    df: pd.DataFrame, id_col: str, suffix: str, value_name: str, sample_ids: list
) -> pd.DataFrame:
    """Generate melted/long format dataframefor a given sample-specific field's prefix (e.g. genotype), for all samples detected"""
    sample_cols = [f"{sample}{suffix}" for sample in sample_ids]
    melted = df.melt(
        id_vars=[id_col],
        value_vars=sample_cols,
        var_name="Sample",
        value_name=value_name,
    )
    melted["Sample"] = melted["Sample"].str.replace(f"{suffix}", "", regex=False)

    return melted

def is_compound_het(
    maternal_haplotype_count, paternal_haplotype_count, unknown_haplotype_count
):
    """Determine if a gene is compound het based on maternal and paternal haplotype counts and unknown haplotype count."""
    if maternal_haplotype_count > 0 and paternal_haplotype_count > 0:
        return "True"
    elif unknown_haplotype_count == 1 and (
        maternal_haplotype_count > 0 or paternal_haplotype_count > 0
    ):
        return "Unknown"
    elif unknown_haplotype_count > 1:
        return "Unknown"
    else:
        return "False"

def count_gene_variants_by_parental_haplotype(
    variant_gt_details: pd.DataFrame, fam_dict: dict
) -> pd.DataFrame:
    """Count variants by parental haplotype for a given gene.

    For each heterozygous variant in the gene, determines whether it was inherited on the maternal
    or paternal haplotype based on parental genotypes. Inheritance is determined by:
    - If mother is heterozygous and father is reference, counts as maternal
    - If father is heterozygous and mother is reference, counts as paternal
    - If one parent genotype is missing, assumes inheritance from heterozygous parent
    - If both parents are reference or missing, counts as unknown

    Args:
        variant_gt_details (pd.DataFrame): DataFrame with variant genotype details including
            zygosity for proband and parents
        fam_dict (dict): Dictionary with family roles

    Returns:
        pd.DataFrame: gene ids to maternal and paternal haplotype counts
    """

    gene_to_haplotype_counts = {}

    for gene in variant_gt_details["Ensembl_gene_id"].unique():
        maternal_haplotype_count = 0
        paternal_haplotype_count = 0
        unknown_haplotype_count = 0
        for variant in variant_gt_details[variant_gt_details["Ensembl_gene_id"] == gene][
            "Variant_id"
        ]:
            proband_zygosity = variant_gt_details[(variant_gt_details["Variant_id"] == variant) & (variant_gt_details["Sample"] == fam_dict["child"])]["Zygosity"].values[0]
            mother_zygosity = variant_gt_details[(variant_gt_details["Variant_id"] == variant) & (variant_gt_details["Sample"] == fam_dict["mother"])]["Zygosity"].values[0]
            father_zygosity = variant_gt_details[(variant_gt_details["Variant_id"] == variant) & (variant_gt_details["Sample"] == fam_dict["father"])]["Zygosity"].values[0]
            if not (
                proband_zygosity == "homozygous_ref"
                or proband_zygosity == "missing"
                or proband_zygosity == "homozygous_alt"
            ):
                if not (
                    mother_zygosity == "homozygous_alt"
                    or father_zygosity == "homozygous_alt"
                ):
                    if not (
                        mother_zygosity == "heterozygous"
                        and father_zygosity == "heterozygous"
                    ):  # if both parents are heterozygous, we are probably not interested in this variant
                        if mother_zygosity == "heterozygous":
                            if (
                                father_zygosity == "homozygous_ref"
                            ):  # clear maternal inheritance
                                maternal_haplotype_count += 1
                            elif (
                                father_zygosity == "missing"
                            ):  # assume maternal inheritance
                                unknown_haplotype_count += 1
                        elif father_zygosity == "heterozygous":
                            if (
                                mother_zygosity == "homozygous_ref"
                            ):  # clear paternal inheritance
                                paternal_haplotype_count += 1
                            elif (
                                mother_zygosity == "missing"
                            ):  # assume paternal inheritance
                                unknown_haplotype_count += 1
                        elif (
                            mother_zygosity == "homozygous_ref"
                            and father_zygosity == "homozygous_ref"
                        ):  # can later use phase set and phase to determine which haplotype to add to?
                            unknown_haplotype_count += 1
                        elif (
                            mother_zygosity == "missing"
                            and father_zygosity == "missing"
                        ):  # can later use phase set and phase to determine which haplotype to add to?
                            unknown_haplotype_count += 1

        gene_to_haplotype_counts[gene] = [
            maternal_haplotype_count,
            paternal_haplotype_count,
            unknown_haplotype_count,
        ]
    gene_to_haplotype_counts = pd.DataFrame.from_dict(
        gene_to_haplotype_counts, orient="index"
    )
    gene_to_haplotype_counts.columns = [
        "maternal_haplotype_count",
        "paternal_haplotype_count",
        "unknown_haplotype_count",
    ]

    return gene_to_haplotype_counts


def SV_comp_het_status(genes, variant_impact, compound_het_status): 
    """
    Determine compound het status for SV variants

    Args:
        genes (list): List of genes impacted by an SV
        compound_het_status (pd.DataFrame): DataFrame of compound het status for genes

    Returns:
        str: Compound het status for the genes impacted by the SV
    """
    gene_compound_het_statuses = []
    if pd.isna(genes): 
        return "."
    elif variant_impact == "intergenic_region":
        return "."
    else:
        genes = re.split(';|-', genes)
        for gene in genes: 
            try:
                gene_compound_het_statuses.append(compound_het_status.loc[gene, "compound_het_status"])
            except:
                gene_compound_het_statuses.append(".")
                
        return ";".join(gene_compound_het_statuses)
    

def SV_SNV_comp_het_status(SV_genes, variant_impact, SV_SNV_dict): 
    """
    Determine compound het status for SV variants
    TODO: check if variants are in cis or in trans using phased GTs

    Args:
        SV_genes (list): List of genes impacted by an SV
        SNV_SV_dict (dict): Dictionary of cross SV SNV compound het status

    Returns:
        str: SNV compound het status for the genes impacted by the SV
    """
    gene_compound_het_statuses = []
    if pd.isna(SV_genes): 
        return "."
    elif variant_impact == "intergenic_region":
        return "."
    else:
        genes = re.split(';|-', SV_genes)
        for gene in genes: 
            try:
                gene_compound_het_statuses.append(SV_SNV_dict[gene])
            except:
                gene_compound_het_statuses.append("False")
                
        return ";".join(gene_compound_het_statuses)


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
    fam_dict = infer_pedigree_roles(args.pedigree)
    proband_id = fam_dict["child"]

    # read in SV file and filter
    SV = pd.read_csv(args.sv, low_memory=False)
    SV["Variant_id"] = SV["CHROM"] + "-" + SV["POS"].astype(str) + "-" + SV["END"].astype(str) +  "-" + SV["SVTYPE"] + "-" + SV["ID"]
    SV_rare_high_impact = SV[(SV["VARIANT"] != "intergenic_region") & (SV["gnomad_maxAF"] <= 0.01)].copy()
    SV_rare_high_impact.rename(columns={"ENSEMBL_GENE": "Ensembl_gene_id"}, inplace=True) # for compatibility with SNV CH functions

    # get per-sample SV genotype status
    # extract all sample IDs based on genotype columns
    sample_ids = sorted(col.rsplit("_", 1)[0] for col in SV.columns if "_GT" in col)
    variant_ps = melt_sample_columns_SV(SV_rare_high_impact, "Variant_id", "_PS", "PS", sample_ids)
    variant_gts = melt_sample_columns_SV(
        SV_rare_high_impact, "Variant_id", "_GT", "GT", sample_ids
    )
    variant_zyg = melt_sample_columns_SV(
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
    compound_het_status_no_parents = determine_compound_het_status_no_parents(
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
    gene_haplotype_counts = count_gene_variants_by_parental_haplotype(
        unknown_variants, fam_dict
    )
    gene_haplotype_counts["compound_het_status"] = gene_haplotype_counts.apply(
        lambda x: is_compound_het(
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

    # TODO:export table of variants in genes with compound het status
    variant_gt_details = variant_gt_details.merge(compound_het_status, on="Ensembl_gene_id", how="left")
    variant_gt_details.to_csv(f"{family}_compound_het_SVs.csv", index=False)

    # add compound het status to SV report
    # TODO: check if variants are in cis or in trans using phased GTs
    # SV["compound_het_status_SV"] = SV_rare_high_impact.apply(lambda x: SV_comp_het_status(x["Ensembl_gene_id"], x["VARIANT"], compound_het_status), axis=1)
    # SNVs = pd.read_csv(args.snvs)
    # SNV_SV_CH = {}
    # for gene in variant_gt_details_proband["Ensembl_gene_id"].unique(): 
    #     if gene in SNVs["Ensembl_gene_id"].values: 
    #         SNV_SV_CH[gene] = "True"
    #     else:
    #         SNV_SV_CH[gene] = "False"
    # SV["compound_het_status_SNV"] = SV.apply(lambda x: SV_SNV_comp_het_status(x["ENSEMBL_GENE"], x["VARIANT"], SNV_SV_CH), axis=1)


main()

 