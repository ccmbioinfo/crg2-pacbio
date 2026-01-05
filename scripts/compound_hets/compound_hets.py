import argparse
from functools import reduce
import pandas as pd
import pyranges as pr
import numpy as np
import re


def parse_promoterAI_score(score):
    """Parse the promoterAI score.
    Args:
        score (str): The promoterAI score.
    Returns:
        float: The parsed promoterAI score.
    """
    if pd.isna(score):
        return None
    if score == "No":
        return None
    try:
        return max([abs(float(i)) for i in score.split(",")])
    except:
        return score


def filter_low_impact_variants(low: pd.DataFrame) -> pd.DataFrame:
    """Select low impact variants that may be damaging

    Args:
        low (pd.DataFrame): DataFrame with low impact variants queried from gemini database

    Returns:
        pd.DataFrame: DataFrame with genic low impact variants that may be damaging (SNVs with SpliceAI >=0.2, CADD >= 10, promoterAI >= 0.1, or indel with no SpliceAI or CADD score)
    """
    # filter out intergenic variants
    low_impact_var_filter = low[
        low["Variation"].isin(
            [
                "intron_variant",
                "3_prime_UTR_variant",
                "5_prime_UTR_variant",
                "non_coding_transcript_exon_variant",
                "mature_miRNA_variant",
                "synonymous_variant",
            ]
        )
    ].copy()
    # parse SpliceAI, CADD, and promoterAI scores
    low_impact_var_filter["SpliceAI_score_parsed"] = low_impact_var_filter[
        "SpliceAI_score"
    ].apply(lambda x: max(x.split("|")[2:6]) if not pd.isna(x) else np.nan)
    low_impact_var_filter["SpliceAI_score_parsed"] = low_impact_var_filter[
        "SpliceAI_score_parsed"
    ].astype(float)
    low_impact_var_filter["Cadd_score"] = low_impact_var_filter["Cadd_score"].astype(
        float
    )
    low_impact_var_filter["promoterAI_score_parsed"] = low_impact_var_filter["promoterAI_score"].apply(parse_promoterAI_score)
    # identify variant type
    low_impact_var_filter["sum_ref_alt_length"] = (
        low_impact_var_filter["Ref"].str.len() + low_impact_var_filter["Alt"].str.len()
    )
    low_impact_var_filter["variant_type"] = low_impact_var_filter[
        "sum_ref_alt_length"
    ].apply(lambda x: "SNV" if x == 2 else "indel")
    # retain SNVs with SpliceAI >=0.2, CADD >= 10, promoterAI >= 0.1, or indel with no SpliceAI or CADD score:
    low_impact_var_filter_scores = low_impact_var_filter[
        (low_impact_var_filter["SpliceAI_score_parsed"] >= 0.2)
        | (low_impact_var_filter["Cadd_score"] >= 10)
        | (low_impact_var_filter["promoterAI_score_parsed"] >= 0.1)
        | (
            (low_impact_var_filter["variant_type"] == "indel")
            & (low_impact_var_filter["SpliceAI_score_parsed"].isna())
            & (low_impact_var_filter["Cadd_score"].isna())
        )  # missing spliceAI and cadd score
    ]

    return low_impact_var_filter_scores


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
    pedigree = pedigree.astype(str)
    # just take ID of first affected child if there are multiple affected children
    try:
        child = pedigree[
            (pedigree["paternal_ID"] != "0") & (pedigree["maternal_ID"] != "0") & (pedigree["phenotype"] == "2")
        ]["individual_ID"].values[0]
    except IndexError:
        # one or both parents are missing, pull first affected individual
        try:
            child =  pedigree[
                (pedigree["phenotype"] == "2")
            ]["individual_ID"].values[0]  
        except IndexError:
            logging.error("No affected individuals found in pedigree")
            raise IndexError("No affected individuals found in pedigree")
    father = pedigree[pedigree["individual_ID"] == child]["paternal_ID"].values[0]
    mother = pedigree[pedigree["individual_ID"] == child]["maternal_ID"].values[0]
    fam_dict = {"child": child, "father": father, "mother": mother}

    return fam_dict


def abstract_gt(ref, alt, gt):
    if "|" not in gt:
        return gt
    else:
        gt_0 = gt.split("|")[0]
        gt_1 = gt.split("|")[1]
        if ref == gt_0:
            return "ref|alt"
        else:
            return "alt|ref"


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
        for variant in variant_gt_details[
            variant_gt_details["Ensembl_gene_id"] == gene
        ]["Variant_id"]:
            proband_zygosity = variant_gt_details[
                (variant_gt_details["Variant_id"] == variant)
                & (variant_gt_details["Sample"] == fam_dict["child"])
            ]["Zygosity"].values[0]
            mother_zygosity = variant_gt_details[
                (variant_gt_details["Variant_id"] == variant)
                & (variant_gt_details["Sample"] == fam_dict["mother"])
            ]["Zygosity"].values[0]
            father_zygosity = variant_gt_details[
                (variant_gt_details["Variant_id"] == variant)
                & (variant_gt_details["Sample"] == fam_dict["father"])
            ]["Zygosity"].values[0]
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
                            elif father_zygosity == "missing":
                                unknown_haplotype_count += 1
                        elif father_zygosity == "heterozygous":
                            if (
                                mother_zygosity == "homozygous_ref"
                            ):  # clear paternal inheritance
                                paternal_haplotype_count += 1
                            elif mother_zygosity == "missing":
                                unknown_haplotype_count += 1
                        elif (
                            mother_zygosity == "homozygous_ref"
                            and father_zygosity == "homozygous_ref"
                        ):
                            unknown_haplotype_count += 1
                        elif (
                            mother_zygosity == "missing"
                            and father_zygosity == "missing"
                        ):
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


def is_compound_het(
    maternal_haplotype_count, paternal_haplotype_count, unknown_haplotype_count
):
    """Determine if a gene is compound het based on maternal and paternal haplotype counts and unknown haplotype count."""
    if maternal_haplotype_count > 0 and paternal_haplotype_count > 0:
        return "TRUE"
    elif unknown_haplotype_count == 1 and (
        maternal_haplotype_count > 0 or paternal_haplotype_count > 0
    ):
        return "UNKNOWN"
    elif unknown_haplotype_count > 1:
        return "UNKNOWN"
    else:
        return "FALSE"


def determine_compound_het_status_no_parents(
    variant_gt_details: pd.DataFrame,
) -> pd.DataFrame:
    """Determine compound het status for a given gene using only long-read phasing information.

    Args:
        variant_to_gene (pd.DataFrame): DataFrame mapping variant IDs to Ensembl gene IDs
        variant_gt_details (pd.DataFrame): DataFrame with variant genotype details

    Returns:
        pd.DataFrame: DataFrame with compound het status for each gene
    """
    gene_to_compound_hets = {}
    for gene in variant_gt_details["Ensembl_gene_id"].unique():
        variants = variant_gt_details[variant_gt_details["Ensembl_gene_id"] == gene]
        num_variants = len(variants)  # total number of het variants in the gene
        true_compound_hets = 0
        unknown_compound_hets = 0
        phase_sets = variants["PS"].unique()
        if num_variants > 1:
            if phase_sets.tolist() == ["."] or phase_sets.tolist() == [
                "TR_overlap"
            ]:  # variants are not phased
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
            compound_het_status = "TRUE"
        else:
            if unknown_compound_hets > 0:
                compound_het_status = "UNKNOWN"
            else:
                compound_het_status = "FALSE"
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


def get_compound_het_variants(
    high_impact: pd.DataFrame,
    variant_to_gene: pd.DataFrame,
    compound_het_status: pd.DataFrame,
    proband_id: str,
    sample_ids: list,
) -> pd.DataFrame:
    # TODO: may want to export all SNVs for cross-variant type CH analysis
    """Get all heterozygous high impactvariants for a given proband in genes annotatedwith compound het status.
    This will be used for cross-variant type CH analysis.

    Args:
        high_impact (pd.DataFrame): DataFrame with high impact variants
        variant_to_gene (pd.DataFrame): DataFrame mapping variant IDs to gene IDs
        compound_het_status (pd.DataFrame): DataFrame with compound het status for each gene
        proband_id (str): ID of the proband
        sample_ids (list): List of sample IDs

    Returns:
        pd.DataFrame: DataFrame with high impact variants that are compound het
    """
    ch_genes = compound_het_status.index.tolist()
    ch_variants = variant_to_gene[variant_to_gene["Ensembl_gene_id"].isin(ch_genes)][
        "Variant_id"
    ].values.tolist()
    ch_variants_df = high_impact[high_impact["Variant_id"].isin(ch_variants)].copy()
    # add zygosity column
    for sample in sample_ids:
        gt_type_col = f"gt_types.{sample}"
        zygosity_col = f"zygosity.{sample}"
        if gt_type_col in ch_variants_df.columns:
            ch_variants_df[zygosity_col] = ch_variants_df[gt_type_col].map(
                lambda x: (
                    "hom"
                    if x == 3
                    else "het" if x == 1 else "-" if x == 0 else "missing"
                )
            )
    ch_variants_df = ch_variants_df[ch_variants_df[f"zygosity.{proband_id}"] == "het"]
    # format dataframe for export
    ch_variants_df["Gnomad_af_popmax"] = ch_variants_df["Gnomad_af_popmax"].replace(
        -1, 0
    )
    ch_variants_df.drop(
        columns=["Variant_id", "variant_type", "sum_ref_alt_length", "GT_qual_max"],
        inplace=True,
    )
    gt_type_drop = [
        col for col in ch_variants_df.columns if col.startswith("gt_types.")
    ]
    gt_phases_drop = [
        col for col in ch_variants_df.columns if col.startswith("gt_phases.")
    ]
    ch_variants_df.drop(columns=gt_type_drop + gt_phases_drop, inplace=True)
    ch_variants_df = ch_variants_df[
        ["Chrom", "Pos", "Ref", "Alt"]
        + ["zygosity." + sample for sample in sample_ids]
        + [
            "Variation",
            "Depth",
            "Quality",
            "Gene",
            "Clinvar",
            "Ensembl_gene_id",
            "Gnomad_af_popmax",
            "Cadd_score",
            "SpliceAI_score",
        ]
        + [f"gts.{sample}" for sample in sample_ids]
        + [f"gt_alt_depths.{sample}" for sample in sample_ids]
        + [f"gt_depths.{sample}" for sample in sample_ids]
        + [f"gt_quals.{sample}" for sample in sample_ids]
        + [f"PS.{sample}" for sample in sample_ids]
    ]

    return ch_variants_df


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


def SV_comp_het_status(genes, variant_impact, compound_het_status):
    """
    Determine compound het status for SV variants

    Args:
        genes (list): List of genes impacted by an SV
        compound_het_status (pd.DataFrame): DataFrame of compound het status for genes

    Returns:
        dict: Dictionary with CH_status and CH_variant_types for the genes impacted by the SV
    """
    gene_compound_het_statuses = []
    gene_compound_het_variant_types = []
    if pd.isna(genes):
        return {"CH_variant_types": ".", "CH_status": "."}
    elif variant_impact == "intergenic_region":
        return {"CH_variant_types": ".", "CH_status": "."}
    else:
        genes = re.split(";|-", genes)
        for gene in genes:
            try:
                gene_compound_het_statuses.append(
                    compound_het_status.loc[gene, "CH_status"]
                )
                gene_compound_het_variant_types.append(
                    compound_het_status.loc[gene, "CH_variant_types"]
                )
            except KeyError:
                gene_compound_het_statuses.append(".")
                gene_compound_het_variant_types.append(".")
        return {
            "CH_variant_types": "|".join(gene_compound_het_variant_types),
            "CH_status": "|".join(gene_compound_het_statuses),
        }  # return a dictionary with CH_status and CH_variant_types


def get_CNV_genotypes(cnv_details):
    if pd.isna(cnv_details):
        return "0/0"
    else:
        return cnv_details.split(":")[3]


def map_gene_symbol_to_ensembl(gene_symbol, ensembl_df):
    """
    Map a gene symbol to an Ensembl gene ID using the provided Ensembl DataFrame.

    Args:
        gene_symbol (str): The gene symbol to map.
        ensembl_df (pd.DataFrame): The Ensembl DataFrame containing gene information.

    Returns:
        str: The Ensembl gene ID if found, otherwise the original gene symbol.
    """
    return ensembl_df[ensembl_df["gene_name"] == gene_symbol]["gene_id"].values[0]


def map_NCBI_ID_to_ensembl(NCBI_ID, ensembl_to_NCBI_df):
    """
    Map an NCBI gene ID to an Ensembl gene ID using the provided Ensembl to NCBI ID dataframe.

    Args:
        NCBI_ID (str): The NCBI gene ID to map.
        ensembl_to_NCBI_df (pd.DataFrame): The Ensembl to NCBI ID dataframe.

    Returns:
        str: The Ensembl gene ID if found, otherwise the gene symbol.
    """
    try:
        return ensembl_to_NCBI_df[ensembl_to_NCBI_df["entrezgene_id"] == int(NCBI_ID)]["ensembl_gene_id"].values[0]
    except:
        return NCBI_ID


def replace_NCBI_IDs_with_Ensembl_IDs(variant_df, ensembl_df, ensembl_to_NCBI_df):
    """
    Replace NCBI gene IDs with Ensembl gene IDs in a variant dataframe.
    Sometimes there are NCBI IDs in the Ensembl_gene_id column. 
    This function replaces the NCBI IDs with the Ensembl IDs to facilitate cross-variant compound het status determination.

    Args:
        variant_df (pd.DataFrame): The variant dataframe to replace NCBI IDs with Ensembl IDs.
        ensembl_df (pd.DataFrame): The Ensembl DataFrame containing gene information.
        ensembl_to_NCBI_df (pd.DataFrame): The Ensembl to NCBI ID dataframe.

    Returns:
        pd.DataFrame: The variant dataframe with NCBI IDs replaced with Ensembl IDs.
    """
    # get genes with NCBI ID instead of Ensembl ID
    no_ensembl_id = variant_df[~variant_df["Ensembl_gene_id"].str.contains("ENSG")]["Gene"].unique()
    ID_dict = {}
    not_mapped = []
    for gene in no_ensembl_id:
        try: 
            ensembl_id = map_gene_symbol_to_ensembl(gene, ensembl_df)
            ID_dict[gene] = ensembl_id
        except:
            not_mapped.append(gene)
    # map NCBI IDs to Ensembl IDs
    for gene in not_mapped:
        NCBI_ID = [id for id in variant_df[variant_df["Gene"] == gene]["Ensembl_gene_id"].values if "ENSG" not in id][0]
        try:
            ensembl_id = map_NCBI_ID_to_ensembl(NCBI_ID, ensembl_to_NCBI_df)
            ID_dict[gene] = ensembl_id
        except:
            not_mapped.append(gene)
        
    # map gene symbols to Ensembl IDs
    mask = ~variant_df["Ensembl_gene_id"].str.contains("ENSG")
    for idx, row in variant_df[mask].iterrows():
        gene = row["Gene"]
        if gene in ID_dict and pd.notnull(ID_dict[gene]):
            variant_df.at[idx, "Ensembl_gene_id"] = ID_dict[gene] 