import argparse
from functools import reduce
import glob
import logging
import os
import pandas as pd
import pyranges as pr
import numpy as np
import re
from typing import Dict, Tuple

from compound_hets import compound_hets
from annotation import annotate

# Constants
ZYGOSITY_HET = "heterozygous"
ZYGOSITY_HOM_ALT = "homozygous_alt"
ZYGOSITY_HOM_REF = "homozygous_ref"
ZYGOSITY_MISSING = "missing"

VARIANT_TYPE_SEQUENCE = "sequence_variant"
VARIANT_TYPE_SV = "SV"
VARIANT_TYPE_CNV = "CNV"

GT_TYPE_HET = 1
GT_TYPE_HOM_ALT = 3
GT_TYPE_MISSING = 2

COMPOUND_HET_STATUS_TRUE = "TRUE"
COMPOUND_HET_STATUS_FALSE = "FALSE"
COMPOUND_HET_STATUS_UNKNOWN = "UNKNOWN"

MISSING_VALUE = "."


def setup_logging() -> logging.Logger:
    """Configure and return logger."""
    logging.basicConfig(
        level=getattr(logging, "INFO", logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    return logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate compound het status across variant types."
    )
    parser.add_argument(
        "--high_med",
        type=str,
        required=True,
        help="Path to HIGH-MED sequence variants TSV",
    )
    parser.add_argument(
        "--low", type=str, required=True, help="Path to LOW sequence variants TSV"
    )
    parser.add_argument(
        "--sequence_variant_report_dir",
        type=str,
        required=True,
        help="Path to directory containing sequence variant report CSV",
    )
    parser.add_argument("--sv", type=str, required=True, help="Path to SV report CSV")
    parser.add_argument("--cnv", type=str, required=True, help="Path to CNV report CSV")
    parser.add_argument(
        "--ensembl", type=str, required=True, help="Path to Ensembl gene CSV"
    )
    parser.add_argument(
        "--ensembl_to_NCBI_df", type=str, required=True, help="Path to Ensembl to NCBI ID CSV"
    )
    parser.add_argument(
        "--pedigree", type=str, required=True, help="Path to pedigree file"
    )
    parser.add_argument("--family", type=str, required=True, help="Family ID")
    return parser.parse_args()


def setup_pedigree(pedigree_path: str, family: str, logger: logging.Logger) -> Tuple[str, str, Dict[str, str]]:
    """Load pedigree and extract family information."""
    fam_dict = compound_hets.infer_pedigree_roles(pedigree_path)
    proband_id = fam_dict["child"]
    logger.info(
        "Loaded pedigree %s (family %s, proband %s)", pedigree_path, family, proband_id
    )
    return family, proband_id, fam_dict


def create_variant_id(chrom: str, pos, ref: str, alt: str) -> str:
    """Create variant ID in the form chr-pos-ref-alt."""
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    return variant_id.replace("chr", "")


def filter_heterozygous_variants(
    df: pd.DataFrame, proband_id: str, fam_dict: Dict[str, str], proband_col: str, 
    mother_col: str, father_col: str, proband_value, hom_alt_value
) -> pd.DataFrame:
    """Filter variants that are het in proband and neither parent is homozygous alternate."""
    if mother_col == "gt_types.0" and father_col == "gt_types.0": # no parents in pedigree
        return df[df[proband_col] == proband_value]
    elif mother_col == "gt_types.0":
        # mother is missing
        return df[
            (df[proband_col] == proband_value)
            & (df[father_col] != hom_alt_value)
        ]
    elif father_col == "gt_types.0":
        # father is missing
        return df[
            (df[proband_col] == proband_value)
            & (df[mother_col] != hom_alt_value)
        ]
    else:
        return df[
            (df[proband_col] == proband_value)
            & (df[mother_col] != hom_alt_value)
            & (df[father_col] != hom_alt_value)
        ]


def extract_sample_ids_from_columns(df: pd.DataFrame, prefix: str) -> list:
    """Extract sample IDs from columns with given prefix."""
    if prefix == "gts.":
        return sorted(
            {
                re.sub(r"^gts\.", "", col)
                for col in df.columns
                if col.startswith("gts.")
            }
        )
    elif prefix.endswith("_GT"):
        return sorted(col.rsplit("_", 1)[0] for col in df.columns if "_GT" in col)
    return []


def merge_melted_dataframes(dfs: list, merge_on: list) -> pd.DataFrame:
    """Successively merge a list of DataFrames on specified columns."""
    return reduce(
        lambda left, right: left.merge(right, on=merge_on, how="outer"),
        dfs,
    )


def map_zygosity_from_gt_type(gt_type: int) -> str:
    """Map GT_type numeric code to zygosity string."""
    if gt_type == 3:
        return ZYGOSITY_HOM_ALT
    elif gt_type == 1:
        return ZYGOSITY_HET
    elif gt_type == 0:
        return ZYGOSITY_HOM_REF
    else:
        return ZYGOSITY_MISSING


def process_sequence_variants(
    high_med_path: str,
    low_path: str,
    proband_id: str,
    fam_dict: Dict[str, str],
    ensembl: pd.DataFrame,
    ensembl_to_NCBI_df: pd.DataFrame,
    logger: logging.Logger,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process sequence variants and return high_impact and variant_gt_details dataframes."""
    logger.info("Reading HIGH-MED sequence variants from %s", high_med_path)
    high_med = pd.read_csv(high_med_path, sep="\t")
    logger.info("Loaded %d HIGH-MED variants", len(high_med))
    
    logger.info("Reading LOW sequence variants from %s", low_path)
    low = pd.read_csv(low_path, sep="\t", low_memory=False)
    logger.info("Loaded %d LOW variants", len(low))

    # Filter low impact variants
    low_impact_var_filter_scores = compound_hets.filter_low_impact_variants(low)
    logger.info(
        "Filtered LOW impact variants to %d potentially high-impact variants",
        len(low_impact_var_filter_scores),
    )
    
    # Combine and deduplicate
    high_impact = pd.concat([low_impact_var_filter_scores, high_med])
    high_impact = high_impact.drop_duplicates(subset=["Variant_id"]).copy()

    # Replace NCBI IDs with Ensembl IDs where neccessary (sequence variant Ensembl_gene_id column may contain NCBI IDs)
    compound_hets.replace_NCBI_IDs_with_Ensembl_IDs(high_impact, ensembl, ensembl_to_NCBI_df)

    # Filter by TG inhouse allele count
    high_impact["TG_LRWGS_ac"].replace(np.nan, 0, inplace=True)
    high_impact["TG_LRWGS_ac"] = high_impact["TG_LRWGS_ac"].astype(int)
    high_impact = high_impact[high_impact["TG_LRWGS_ac"] < 20]
    
    # Create variant ID
    high_impact["Variant_id"] = high_impact.apply(
        lambda x: create_variant_id(x["Chrom"], x["Pos"], x["Ref"], x["Alt"]), axis=1
    )

    # Filter heterozygous variants
    high_impact = filter_heterozygous_variants(
        high_impact,
        proband_id,
        fam_dict,
        f"gt_types.{proband_id}",
        f"gt_types.{fam_dict['mother']}",
        f"gt_types.{fam_dict['father']}",
        GT_TYPE_HET,
        GT_TYPE_HOM_ALT,
    )
    logger.info(
        "Retained %d heterozygous high-impact variants for proband %s",
        len(high_impact),
        proband_id,
    )

    # Extract sample IDs
    sample_ids = extract_sample_ids_from_columns(high_impact, "gts.")
    print(sample_ids)
    
    # Split PS column
    if len(sample_ids) > 1:
        ps_split_cols = high_impact["PS"].str.split(",", expand=True)
        for idx, sample in enumerate(sample_ids):
            high_impact[f"PS.{sample}"] = ps_split_cols[idx]
    else:
        high_impact[f"PS.{sample_ids[0]}"] = high_impact["PS"]
    high_impact.drop(columns=["PS"], inplace=True)

    logger.debug("Extracted %d sample IDs from sequence variant table", len(sample_ids))

    # Melt sample columns
    variant_ps = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "PS.", "PS", sample_ids
    )
    variant_gts = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gts.", "GT", sample_ids
    )
    variant_gt_type = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gt_types.", "GT_type", sample_ids
    )
    variant_qual = compound_hets.melt_sample_columns(
        high_impact, "Variant_id", "gt_quals.", "GT_qual", sample_ids
    )

    # Merge melted dataframes
    dfs = [variant_ps, variant_gts, variant_gt_type, variant_qual]
    sequence_variant_gt_details = merge_melted_dataframes(
        dfs, ["Variant_id", "Ref", "Alt", "Sample"]
    )
    
    # Fix missing genotypes
    sequence_variant_gt_details.loc[
        sequence_variant_gt_details["GT"] == "./.", "GT_type"
    ] = GT_TYPE_MISSING
    
    # Map zygosity
    sequence_variant_gt_details["Zygosity"] = sequence_variant_gt_details[
        "GT_type"
    ].map(map_zygosity_from_gt_type)
    
    # Add gene information
    sequence_variant_to_gene = high_impact[
        ["Variant_id", "Ensembl_gene_id"]
    ].drop_duplicates()
    sequence_variant_gt_details = sequence_variant_gt_details.merge(
        sequence_variant_to_gene, on="Variant_id", how="left"
    )
    logger.debug(
        "Sequence variant details dataframe contains %d rows",
        len(sequence_variant_gt_details),
    )

    # Abstract genotype
    sequence_variant_gt_details["GT_abstracted"] = sequence_variant_gt_details.apply(
        lambda x: compound_hets.abstract_gt(x["Ref"], x["Alt"], x["GT"]), axis=1
    )

    return high_impact, sequence_variant_gt_details


def process_structural_variants(
    sv_path: str,
    proband_id: str,
    fam_dict: Dict[str, str],
    variant_type: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """Process structural variants and return variant_gt_details dataframe."""
    logger.info("Reading %s report from %s", variant_type, sv_path)
    SV = pd.read_csv(sv_path, low_memory=False)
    
    # Filter rare genic variants
    SV_rare_high_impact = SV[
        (SV["VARIANT"] != "intergenic_region") & (SV["gnomad_maxAF"] <= 0.01) & (SV["TG_nhomalt_max"] <= 5)
    ].copy()
    logger.info(
        "Identified %d rare genic %s (gnomAD AF <= 0.01)",
        variant_type,
        len(SV_rare_high_impact),
    )
    
    # Create variant ID
    if variant_type == VARIANT_TYPE_SV: # SVs
        SV_rare_high_impact["Variant_id"] = (
            SV_rare_high_impact["CHROM"].astype(str)
            + "-"
            + SV_rare_high_impact["POS"].astype(str)
            + "-"
            + SV_rare_high_impact["END"].astype(str)
            + "-"
            + SV_rare_high_impact["SVTYPE"]
            + "-"
            + SV_rare_high_impact["ID"]
        )
    else: # CNVs
        SV_rare_high_impact["Variant_id"] = (
            SV_rare_high_impact["CHROM"].astype(str)
            + "-"
            + SV_rare_high_impact["POS"].astype(str)
            + "-"
            + SV_rare_high_impact["END"].astype(str)
            + "-"
            + SV_rare_high_impact["SVTYPE"]
        )
    SV_rare_high_impact.rename(
        columns={"ENSEMBL_GENE": "Ensembl_gene_id"}, inplace=True
    )
    
    # Filter heterozygous variants
    if fam_dict.get("mother") == "0" and fam_dict.get("father") == "0": # no parents in pedigree
        SV_rare_high_impact = SV_rare_high_impact[SV_rare_high_impact[f"{proband_id}_zyg"] == "het"]
    elif fam_dict.get("mother") == "0":
        # mother is missing
        SV_rare_high_impact = SV_rare_high_impact[
            (SV_rare_high_impact[f"{proband_id}_zyg"] == "het")
            & (SV_rare_high_impact[f"{fam_dict['father']}_zyg"] != "hom")
        ]
    elif fam_dict.get("father") == "0":
        # father is missing
        SV_rare_high_impact = SV_rare_high_impact[
            (SV_rare_high_impact[f"{proband_id}_zyg"] == "het")
            & (SV_rare_high_impact[f"{fam_dict['mother']}_zyg"] != "hom")
        ]
    else:
        SV_rare_high_impact = SV_rare_high_impact[
            (SV_rare_high_impact[f"{proband_id}_zyg"] == "het")
            & (SV_rare_high_impact[f"{fam_dict['mother']}_zyg"] != "hom")
            & (SV_rare_high_impact[f"{fam_dict['father']}_zyg"] != "hom")
        ]
    logger.info(
        "Retained %d heterozygous %s for proband %s",
        variant_type,
        len(SV_rare_high_impact),
        proband_id,
    )

    # Extract sample IDs
    sample_ids = extract_sample_ids_from_columns(SV, "_GT")
    
    # Melt sample columns
    if variant_type == VARIANT_TYPE_CNV:
        # CNVs don't have PS columns, create dummy dataframe with PS values as "."
        variant_ids = SV_rare_high_impact["Variant_id"].unique()
        variant_ps = pd.DataFrame({
            "Variant_id": np.repeat(variant_ids, len(sample_ids)),
            "Sample": np.tile(sample_ids, len(variant_ids)),
            "PS": MISSING_VALUE
        })
    else:
        variant_ps = compound_hets.melt_sample_columns_SV(
            SV_rare_high_impact, "Variant_id", "_PS", "PS", sample_ids
        )
    variant_gts = compound_hets.melt_sample_columns_SV(
        SV_rare_high_impact, "Variant_id", "_GT", "GT", sample_ids
    )
    variant_zyg = compound_hets.melt_sample_columns_SV(
        SV_rare_high_impact, "Variant_id", "_zyg", "Zygosity", sample_ids
    )
    
    # Merge melted dataframes
    dfs = [variant_ps, variant_gts, variant_zyg]
    SV_gt_details = merge_melted_dataframes(dfs, ["Variant_id", "Sample"])
    
    # Process gene information
    SV_to_gene = SV_rare_high_impact[
        ["Variant_id", "Ensembl_gene_id"]
    ].drop_duplicates()
    SV_to_gene = SV_to_gene[SV_to_gene["Ensembl_gene_id"].notna()]
    
    # Split and explode gene IDs
    SV_to_gene["Ensembl_gene_id"] = SV_to_gene["Ensembl_gene_id"].str.split(";")
    SV_to_gene = SV_to_gene.explode("Ensembl_gene_id")
    SV_to_gene["Ensembl_gene_id"] = SV_to_gene["Ensembl_gene_id"].str.split("-")
    SV_to_gene = SV_to_gene.explode("Ensembl_gene_id")
    
    # Add gene information
    SV_gt_details = SV_gt_details.merge(SV_to_gene, on="Variant_id", how="left")
    logger.debug("SV details dataframe contains %d rows", len(SV_gt_details))
    
    # Recode zygosity
    SV_gt_details.replace(
        {
            "het": ZYGOSITY_HET,
            "hom": ZYGOSITY_HOM_ALT,
            "-": ZYGOSITY_HOM_REF,
            "./.": ZYGOSITY_MISSING,
        },
        inplace=True,
    )
    SV_gt_details["GT_abstracted"] = SV_gt_details["GT"].replace(
        {"0|1": "ref|alt", "1|0": "alt|ref"}
    )

    return SV_gt_details


def process_cnvs(
    cnv_path: str,
    ensembl_path: str,
    proband_id: str,
    fam_dict: Dict[str, str],
    sample_ids: list,
    logger: logging.Logger,
) -> pd.DataFrame:
    """Process CNVs and return variant_gt_details dataframe."""
    logger.info("Reading CNV report from %s", cnv_path)
    CNV = pd.read_csv(cnv_path, low_memory=False)
    
    # Create variant ID
    CNV["Variant_id"] = (
        CNV["CHROM"].astype(str)
        + "-"
        + CNV["START"].astype(str)
        + "-"
        + CNV["END"].astype(str)
        + "-"
        + CNV["SVTYPE"].astype(str)
    )
    
    # Filter rare variants
    CNV_rare = CNV[
        (CNV["pacBioPctFreq_50pctRecOvlp"] < 5) & (CNV[f"gene_symbol"] != ".")
    ].copy()
    
    # Get genotypes
    for sample in sample_ids:
        CNV_rare[f"{sample}_genotype"] = CNV_rare[f"{sample}_SV_DETAILS"].apply(
            compound_hets.get_CNV_genotypes
        )
    
    # Filter heterozygous variants
    if fam_dict.get("mother") == "0" and fam_dict.get("father") == "0": # no parents in pedigree
        CNV_rare = CNV_rare[CNV_rare[f"{proband_id}_genotype"] == "0/1"]
    elif fam_dict.get("mother") == "0":
        # mother is missing
        CNV_rare = CNV_rare[
            (CNV_rare[f"{proband_id}_genotype"] == "0/1")
            & (CNV_rare[f"{fam_dict['father']}_genotype"] != "1/1")
        ]
    elif fam_dict.get("father") == "0":
        # father is missing
        CNV_rare = CNV_rare[
            (CNV_rare[f"{proband_id}_genotype"] == "0/1")
            & (CNV_rare[f"{fam_dict['mother']}_genotype"] != "1/1")
        ]
    else:
        CNV_rare = CNV_rare[
            (CNV_rare[f"{proband_id}_genotype"] == "0/1")
            & (CNV_rare[f"{fam_dict['mother']}_genotype"] != "1/1")
            & (CNV_rare[f"{fam_dict['father']}_genotype"] != "1/1")
        ]

    # Annotate with Ensembl genes
    ensembl = pd.read_csv(ensembl_path, low_memory=False)
    ensembl["Chromosome"] = "chr" + ensembl["Chromosome"].astype(str)
    ensembl = pr.PyRanges(ensembl)
    CNV_pr = pr.PyRanges(
        CNV_rare[["CHROM", "START", "END", "Variant_id"]].rename(
            columns={"CHROM": "Chromosome", "START": "Start", "END": "End"}
        )
    )
    CNV_ensembl = annotate.annotate_genes(CNV_pr, ensembl)
    CNV_ensembl = CNV_ensembl.rename(columns={"gene_id": "Ensembl_gene_id"})
    CNV_ensembl = CNV_ensembl[["Variant_id", "Ensembl_gene_id"]].drop_duplicates()

    # Melt sample-level genotype columns
    gt_cols = [c for c in CNV_rare.columns if "genotype" in c]
    CNV_rare = CNV_rare[["Variant_id"] + gt_cols]
    CNV_gt_details = compound_hets.melt_sample_columns_SV(
        CNV_rare, "Variant_id", "_genotype", "Genotype", sample_ids
    )
    CNV_gt_details["PS"] = MISSING_VALUE  # CNVs are not phased
    CNV_gt_details["GT_abstracted"] = CNV_gt_details["Genotype"]
    CNV_gt_details["Zygosity"] = CNV_gt_details["Genotype"].apply(
        lambda x: (
            ZYGOSITY_HET
            if x == "0/1"
            else ZYGOSITY_HOM_ALT if x == "1/1" else ZYGOSITY_HOM_REF
        )
    )
    
    # Add gene information
    CNV_gt_details = CNV_gt_details.merge(CNV_ensembl, on="Variant_id", how="left")
    CNV_gt_details = CNV_gt_details[~CNV_gt_details["Ensembl_gene_id"].isna()]
    CNV_gt_details["Variant_type"] = VARIANT_TYPE_CNV

    return CNV_gt_details


def determine_compound_het_status(
    all_variants: pd.DataFrame,
    proband_id: str,
    fam_dict: Dict[str, str],
    logger: logging.Logger,
) -> pd.DataFrame:
    """Determine compound het status for all variants."""
    # Get variants that are heterozygous in the proband
    all_variants_proband = all_variants[
        (all_variants["Sample"] == proband_id)
        & (all_variants["Zygosity"] == ZYGOSITY_HET)
    ].copy()
    logger.info(
        "Combined variant table contains %d rows; %d heterozygous proband variants",
        len(all_variants),
        len(all_variants_proband),
    )
    
    # Phase variants using proband phasing only
    logger.info("Attempting to phase variants using proband phasing only")
    compound_het_status_no_parents = (
        compound_hets.determine_compound_het_status_no_parents(all_variants_proband)
    )
    
    # Check if parents exist in pedigree before attempting parental phasing
    if fam_dict.get("mother") == "0" or fam_dict.get("father") == "0":
        logger.info(
            "Skipping parental phasing: proband is missing one or both parents"
        )
        return compound_het_status_no_parents
    
    # Attempt parental phasing for unknown genes
    unknown_gene = compound_het_status_no_parents[
        compound_het_status_no_parents["compound_het_status"] == COMPOUND_HET_STATUS_UNKNOWN
    ]
    unknown_variants = (
        all_variants[
            all_variants["Ensembl_gene_id"].isin(unknown_gene["Ensembl_gene_id"])
        ]
        .copy()
        .drop_duplicates()
    )
    
    logger.info(
        "Attempting parental phasing for %d genes with unknown status",
        len(unknown_gene),
    )
    logger.debug(
        "Unknown genes: %s", ", ".join(sorted(unknown_gene["Ensembl_gene_id"].unique()))
    )
    
    compound_het_status_parents = (
        compound_hets.count_gene_variants_by_parental_haplotype(
            unknown_variants, fam_dict
        )
    )
    compound_het_status_parents["compound_het_status"] = (
        compound_het_status_parents.apply(
            lambda x: compound_hets.is_compound_het(
                x["maternal_haplotype_count"],
                x["paternal_haplotype_count"],
                x["unknown_haplotype_count"],
            ),
            axis=1,
        )
    )
    
    # Replace unknown status with parental phasing results
    for gene in compound_het_status_parents.index:
        compound_het_status_no_parents.loc[
            compound_het_status_no_parents["Ensembl_gene_id"] == gene,
            "compound_het_status",
        ] = compound_het_status_parents.loc[gene, "compound_het_status"]
    
    return compound_het_status_no_parents


def annotate_reports(
    all_variants: pd.DataFrame,
    sequence_variant_report_dir: str,
    sv_path: str,
    cnv_path: str,
    family: str,
    ensembl_to_NCBI_df: pd.DataFrame,
    logger: logging.Logger,
) -> None:
    """Annotate variant reports with compound het status."""
    # Create gene-level compound het status table
    gene_CH_status = all_variants[
        ["Ensembl_gene_id", "Variant_type", "compound_het_status"]
    ].drop_duplicates()
    gene_CH_status = (
        gene_CH_status.groupby("Ensembl_gene_id")
        .agg({"Variant_type": ";".join, "compound_het_status": "first"})
        .reset_index()
    )
    gene_CH_status.rename(
        columns={
            "Variant_type": "CH_variant_types",
            "compound_het_status": "CH_status",
        },
        inplace=True,
    )

    # Annotate sequence variant report
    sequence_variant_report_path = glob.glob(f"{sequence_variant_report_dir}/*.wgs.coding.*.csv")[0]
    sequence_variant_report = pd.read_csv(sequence_variant_report_path)
    sequence_variant_report = sequence_variant_report.merge(
        gene_CH_status, on="Ensembl_gene_id", how="left"
    )
    sequence_variant_report = sequence_variant_report.fillna(MISSING_VALUE)
    sequence_variant_report.to_csv(
        f"reports/{family}.wgs.coding.CH.csv", index=False
    )
    logger.info(
        "Wrote annotated sequence variant report to %s.wgs.coding.CH.csv",
        family,
    )

    # Annotate SV report
    SV = pd.read_csv(sv_path, low_memory=False)
    CH_status_series = SV.apply(
        lambda x: compound_hets.SV_comp_het_status(
            x["ENSEMBL_GENE"], x["VARIANT"], gene_CH_status.set_index("Ensembl_gene_id"), ensembl_to_NCBI_df
        ),
        axis=1,
    )

    SV["CH_variant_types"] = [
        (
            d["CH_variant_types"]
            if isinstance(d, dict) and "CH_variant_types" in d
            else MISSING_VALUE
        )
        for d in CH_status_series
    ]
    SV["CH_status"] = [
        d["CH_status"] if isinstance(d, dict) and "CH_status" in d else MISSING_VALUE
        for d in CH_status_series
    ]
    SV = SV.fillna(MISSING_VALUE)
    SV.to_csv(f"reports/{family}.sv.CH.csv", index=False)
    logger.info("Wrote annotated SV report to %s.pbsv.CH.csv", family)

    # Annotate CNV report
    CNV = pd.read_csv(cnv_path, low_memory=False)
    CH_status_series = CNV.apply(
        lambda x: compound_hets.SV_comp_het_status(
            x["ENSEMBL_GENE"], x["VARIANT"], gene_CH_status.set_index("Ensembl_gene_id"), ensembl_to_NCBI_df
        ),
        axis=1, 
    )

    CNV["CH_variant_types"] = [
        (
            d["CH_variant_types"]
            if isinstance(d, dict) and "CH_variant_types" in d
            else MISSING_VALUE
        )
        for d in CH_status_series
    ]
    CNV["CH_status"] = [
        d["CH_status"] if isinstance(d, dict) and "CH_status" in d else MISSING_VALUE
        for d in CH_status_series
    ]
    CNV = CNV.fillna(MISSING_VALUE)
    CNV.to_csv(f"reports/{family}.cnv.CH.csv", index=False)
    logger.info("Wrote annotated CNV report to %s.cnv.CH.csv", family)


def main():
    """Main function to orchestrate compound het annotation pipeline."""
    args = parse_arguments()
    logger = setup_logging()
    
    logger.info("Starting compound het annotation")
    if not os.path.exists("reports"):
        os.makedirs("reports")
    family, proband_id, fam_dict = setup_pedigree(args.pedigree, args.family, logger)
    
    # Process sequence variants
    ensembl = pd.read_csv(args.ensembl, low_memory=False)
    ensembl_to_NCBI_df = pd.read_csv(args.ensembl_to_NCBI_df, low_memory=False)
    high_impact, sequence_variant_gt_details = process_sequence_variants(
        args.high_med, args.low, proband_id, fam_dict, ensembl, ensembl_to_NCBI_df, logger
    )
    
    # Extract sample IDs from sequence variants
    sample_ids = extract_sample_ids_from_columns(high_impact, "gts.")
    
    # Process structural variants
    SV_gt_details = process_structural_variants(
        args.sv, proband_id, fam_dict, VARIANT_TYPE_SV, logger
    )
    
    # Process CNVs
    CNV_gt_details = process_structural_variants(
        args.cnv, proband_id, fam_dict, VARIANT_TYPE_CNV, logger
    )
    
    # Combine all variant types
    sequence_variant_gt_details = sequence_variant_gt_details[
        ["Variant_id", "Sample", "PS", "GT_abstracted", "Zygosity", "Ensembl_gene_id"]
    ]
    sequence_variant_gt_details["Variant_type"] = VARIANT_TYPE_SEQUENCE
    
    SV_gt_details = SV_gt_details[
        ["Variant_id", "Sample", "PS", "GT_abstracted", "Zygosity", "Ensembl_gene_id"]
    ]
    SV_gt_details["Variant_type"] = VARIANT_TYPE_SV
    
    CNV_gt_details = CNV_gt_details[
        ["Variant_id", "Sample", "PS", "GT_abstracted", "Zygosity", "Ensembl_gene_id"]
    ]
    CNV_gt_details["Variant_type"] = VARIANT_TYPE_CNV
    
    all_variants = pd.concat(
        [sequence_variant_gt_details, SV_gt_details, CNV_gt_details]
    )
    all_variants = all_variants[~all_variants["Ensembl_gene_id"].isna()]
    
    # Determine compound het status
    compound_het_status = determine_compound_het_status(
        all_variants, proband_id, fam_dict, logger
    )
    
    # Add compound het status to variants
    all_variants = all_variants.merge(
        compound_het_status, on="Ensembl_gene_id", how="left"
    ).sort_values(by=["Ensembl_gene_id", "Sample"])
    
    # Create wide variant table
    ensembl["Chromosome"] = "chr" + ensembl["Chromosome"].astype(str)
    ensembl_pr = pr.PyRanges(ensembl)
    
    wide_variants = all_variants.pivot_table(
        index=["Variant_id", "Variant_type", "Ensembl_gene_id", "compound_het_status"],
        columns="Sample",
        values=["PS", "GT_abstracted", "Zygosity"],
        aggfunc="first",
    )
    wide_variants.columns = [
        f"{field}_{sample}" for field, sample in wide_variants.columns
    ]
    wide_variants = wide_variants.reset_index()
    
    # Add sequence variant annotations
    high_impact_subset = high_impact[
        [
            "Variant_id",
            "Variation",
            "Clinvar",
            "Gnomad_af_grpmax",
            "Cadd_score",
            "SpliceAI_score",
            "promoterAI_score",
            "Nucleotide_change_ensembl",
        ]
    ]
    high_impact_subset["Gnomad_af_grpmax"] = high_impact_subset["Gnomad_af_grpmax"].replace(-1, 0)
    wide_variants = wide_variants.merge(high_impact_subset, on="Variant_id", how="left")
    
    # Add gene name
    ensembl_df = ensembl_pr.df.drop_duplicates(subset="gene_id")
    wide_variants["Gene_name"] = wide_variants["Ensembl_gene_id"].map(
        ensembl_df.set_index("gene_id")["gene_name"]
    )
    
    # Format for export
    zygosity_cols = [c for c in wide_variants.columns if "Zygosity" in c]
    PS_cols = [c for c in wide_variants.columns if "PS" in c]
    GT_cols = [c for c in wide_variants.columns if "GT_abstracted" in c]
    wide_variants.rename(columns={"compound_het_status": "CH_status"}, inplace=True)
    wide_variants = wide_variants[
        [
            "Variant_id",
            "Variant_type",
            "Gene_name",
            "Ensembl_gene_id",
            "CH_status",
        ]
        + zygosity_cols
        + [
            "Variation",
            "Clinvar",
            "Gnomad_af_grpmax",
            "Cadd_score",
            "SpliceAI_score",
            "promoterAI_score",
            "Nucleotide_change_ensembl",
        ]
        + GT_cols
        + PS_cols
    ].sort_values(by=["Ensembl_gene_id"])
    
    wide_variants = wide_variants.replace(
        {
            ZYGOSITY_HET: "het",
            ZYGOSITY_HOM_ALT: "hom",
            ZYGOSITY_HOM_REF: "-",
            np.nan: MISSING_VALUE,
        }
    )
    
    wide_variants.to_csv(f"reports/{family}.compound.het.status.csv", index=False)
    logger.info(
        "Wrote combined variant annotations to %s.compound.het.status.csv", family
    )
    
    # Annotate reports
    annotate_reports(
        all_variants,
        args.sequence_variant_report_dir,
        args.sv,
        args.cnv,
        family,
        ensembl_to_NCBI_df,
        logger,
    )
    
    logger.info("Compound het annotation complete")


if __name__ == "__main__":
    main()
