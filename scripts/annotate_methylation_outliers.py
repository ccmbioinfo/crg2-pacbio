import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr
import sys
from typing import Optional

from annotation.annotate import (
    annotate_genes,
    gene_set,
    add_constraint,
    prepare_OMIM,
    annotate_OMIM,
    add_hpo,
    group_by_gene,
    annotate_reg_regions,
    group_by_greendb,
)

def find_closest_exon(outliers: pd.DataFrame, exons: pr.PyRanges) -> pd.DataFrame:
    """Find the closest exon to an outlier region"""
    outliers_pr = pr.PyRanges(outliers[["Chromosome", "Start", "End"]])
    outliers_exon = outliers_pr.nearest(
        exons, strandedness=False, apply_strand_suffix=False, suffix="_exon"
    ).df
    outliers_exon["closest_exon_gene"] = outliers_exon.apply(
        lambda x: x["gene_name"] if not pd.isna(x["gene_name"]) else x["gene_id"],
        axis=1,
    )
    outliers_exon.drop(
        columns=[
            "Start_exon",
            "End_exon",
            "Strand",
            "gene_name",
            "gene_id",
            "gene_biotype",
            "Feature",
            "Distance",
        ],
        inplace=True,
    )

    return outliers_exon


def SVs_to_pr(svs: str, sample: str) -> pr.PyRanges:
    """
    Filter SVs for sample of interest and convert to PyRanges object
    """
    svs = pd.read_csv(svs, low_memory=False)
    svs.rename(
        columns={"CHROM": "Chromosome", "POS": "Start", "END": "End"}, inplace=True
    )
    # filter out variants that are hom ref in this sample
    svs = svs[svs[f"{sample}_GT"] != "0/0"]
    svs = svs[["Chromosome", "Start", "End", "SVTYPE", f"{sample}_GT", "cmh_maxAF"]]
    svs_pr = pr.PyRanges(svs)

    return svs_pr


def CNVs_to_pr(cnvs: str, sample: str) -> pr.PyRanges:
    """
    Filter CNVs for sample of interest and convert to PyRanges object
    """
    cnvs = pd.read_csv(cnvs, sep="\t")
    cnvs.rename(
        columns={"CHROM": "Chromosome", "START": "Start", "END": "End"}, inplace=True
    )
    cnvs = cnvs[["Chromosome", "Start", "End", "SVTYPE", f"{sample}|GT"]]
    cnvs["Chromosome"] = cnvs["Chromosome"].astype(str)
    cnvs["Chromosome"] = cnvs["Chromosome"].str.replace("chr", "")
    cnvs_pr = pr.PyRanges(cnvs)

    return cnvs_pr


def TRs_to_pr(trs: str, sample: str) -> pr.PyRanges:
    """
    Filter Tandem repeat for sample of interest and z-score >= 3 and convert to PyRanges object
    """
    trs = pd.read_csv(trs)
    trs = trs[(trs["sample"] == sample) & (trs["z_score_len"] >= 3)]
    trs["Chromosome"] = trs["case_trid"].str.split("_").str[0].str.replace("chr", "")
    trs["Start"] = trs["case_trid"].str.split("_").str[1].astype(int)
    trs["End"] = trs["case_trid"].str.split("_").str[2].astype(int)
    trs["motif"] = trs["case_trid"].str.split("_").str[3]
    trs = trs[["Chromosome", "Start", "End", "z_score_len", "motif"]]
    trs_pr = pr.PyRanges(trs)

    return trs_pr


def small_variants_to_pr(small_variants: str) -> pr.PyRanges:
    """
    Convert pre-filtered small variants to PyRanges object
    """
    small_variants = pd.read_csv(
        small_variants,
        header=None,
        sep="\t",
        names=["Chromosome", "Start", "End", "Ref", "Alt", "GT", "gnomAD_AF_popmax"],
    )
    small_variants["Chromosome"] = small_variants["Chromosome"].astype(str)
    small_variants["Chromosome"] = small_variants["Chromosome"].str.replace("chr", "")
    small_variants_pr = pr.PyRanges(small_variants)

    return small_variants_pr


def find_closest_variant(
    outliers_pr: pr.PyRanges, variants_pr: pr.PyRanges, sample: str, suffix: str
) -> pd.DataFrame:
    """
    Find the closest variant to each methylation outlier using PyRanges nearest function
    """
    outliers_pr_variant = outliers_pr.nearest(
        variants_pr, strandedness=False, suffix=f"_{suffix}"
    )
    outliers_variant = outliers_pr_variant.df
    outliers_variant[f"closest_{suffix}"] = outliers_variant.apply(
        lambda x: collapse_variant_anno(x, suffix, sample), axis=1
    )
    outliers_variant.drop(columns=[f"Start_{suffix}", f"End_{suffix}"], inplace=True)
    if suffix == "SV":
        outliers_variant.drop(
            columns=["SVTYPE", f"{sample}_GT", "cmh_maxAF", "Distance"], inplace=True
        )
    elif suffix == "CNV":
        outliers_variant.drop(
            columns=["SVTYPE", f"{sample}|GT", "Distance"], inplace=True
        )
    elif suffix == "TR":
        outliers_variant.drop(
            columns=["z_score_len", "motif", "Distance"], inplace=True
        )
    elif suffix == "SMV":
        outliers_variant.drop(
            columns=["Ref", "Alt", "GT", "gnomAD_AF_popmax", "Distance"], inplace=True
        )

    return outliers_variant


def collapse_variant_anno(row: pd.Series, variant: str, sample: str) -> str:
    # TODO: need to figure out filtering for SV rarity
    # CNVs: no filter?
    # TRs: Z-score >=3
    # SMVs: gnomAD popmax AF < 0.01
    if row.Distance > 1000 and variant != "SMV":
        return np.nan
    elif row.Distance > 50 and variant == "SMV":
        return np.nan
    else:
        if variant == "SV":
            GT_col = [i for i in row.index if "GT" in i and sample in i][0]
            variant = (
                "SV"
                + ":"
                + row.Chromosome
                + ":"
                + str(row.Start_SV)
                + "-"
                + str(row.End_SV)
                + "-"
                + row.SVTYPE
                + "-"
                + row[GT_col]
                + "-AF:"
                + str(row.cmh_maxAF)
                + "-dist:"
                + str(row.Distance)
            )
        elif variant == "CNV":
            GT_col = [i for i in row.index if "GT" in i and sample in i][0]
            variant = (
                "CNV"
                + ":"
                + row.Chromosome
                + ":"
                + str(row.Start_CNV)
                + "-"
                + str(row.End_CNV)
                + "-"
                + row.SVTYPE
                + "-"
                + row[GT_col]
                + "-"
                + "dist:"
                + str(row.Distance)
            )
        elif variant == "TR":
            variant = (
                "TR"
                + ":"
                + row.Chromosome
                + ":"
                + str(row.Start_TR)
                + "-"
                + str(row.End_TR)
                + "-"
                + row.motif
                + "-"
                + "TR"
                + "-"
                + "zscore:"
                + str(row.z_score_len)
                + "-"
                + "dist:"
                + str(row.Distance)
            )
        elif variant == "SMV":
            variant = (
                "small_var"
                + ":"
                + row.Chromosome
                + ":"
                + str(row.Start_SMV)
                + "-"
                + str(row.End_SMV)
                + "-"
                + row.Ref
                + "-"
                + row.Alt
                + "-"
                + row.GT
                + "-"
                + "AF:"
                + str(row.gnomAD_AF_popmax)
                + "-"
                + "dist:"
                + str(row.Distance)
            )
        return variant


def merge_adjacent_outliers(outliers: pd.DataFrame) -> pd.DataFrame:
    """
    Merge adjacent outliers into a single row
    """
    outliers_minimal = pr.PyRanges(outliers[["Chromosome", "Start", "End"]])
    merged_outliers = outliers_minimal.cluster(
        count=True
    )  # use pyranges cluster function to merge adjacent outliers
    outliers_merged_intervals = outliers.merge(
        merged_outliers.df, on=["Chromosome", "Start", "End"], how="left"
    ).sort_values(
        by="Count", ascending=False
    )  # merge outliers with cluster counts
    outliers_merged_grouped = outliers_merged_intervals.groupby(["Cluster"]).agg(
        {  # aggregate by cluster
            "Count": "first",  # count of adjacent outliers
            "Chromosome": "first",
            "Start": "min",
            "End": "max",
            "summary_label": ";".join,
            "compare_label": ";".join,
            "category_pop_count": "max",
            "category_pop_freq": "max",
            "asm_fishers_pvalue": "min",
            "mean_hap1_methyl": lambda x: round(x.mean(), 2),
            "mean_hap2_methyl": lambda x: round(x.mean(), 2),
            "mean_meth_delta": "max",
            "mean_abs_meth_delta_zscore": lambda x: (
                min(x) if min(x) < 0 else max(x)
            ),  # if min zscore is negative, use min, otherwise use max
            "mean_combined_methyl": lambda x: round(x.mean(), 2),
            "mean_combined_methyl_zscore": lambda x: min(x) if min(x) < 0 else max(x),
            "num_phased_cpgs": "sum",
            "num_unphased_cpgs": "sum",
            "num_partial_cpgs": "sum",
            "mean_coverage": lambda x: round(x.mean(), 2),
        }
    )

    return outliers_merged_grouped


def main(
    outliers: str,
    out_file: str,
    ensembl: str,
    constraint: str,
    omim: str,
    SV_path: str,
    CNV_path: str,
    TR_outlier_path: str,
    SMV_path: str,
    coverage: Optional[str] = None,
    hpo: Optional[str] = None,
) -> None:

    print("Loading and filtering methylation outliers")
    sample = outliers.split("/")[-1].split(".")[0]
    outliers = pd.read_csv(outliers, sep="\t")
    outliers.rename(
        {"chrom": "Chromosome", "start": "Start", "end": "End"}, axis=1, inplace=True
    )
    try:
        coverage = pd.read_csv(
            coverage,
            compression="gzip",
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "cpg_ID", "mean_coverage"],
        )
        # annotate regions with average coverage from mosdepth
        outliers = outliers.merge(
            coverage,
            right_on=["chrom", "start", "end"],
            left_on=["Chromosome", "Start", "End"],
            how="left",
        )
    except:
        print("No coverage file provided")

    outliers["Chromosome"] = outliers["Chromosome"].str.replace("chr", "")

    try:
        outliers = outliers[
            (outliers["compare_label"] != "Uncategorized")
            & (outliers["compare_label"] != "InsufficientData")
        ]
    except KeyError:  # cohort comparison
        outliers = outliers[
            (outliers["summary_comparison"] != "Uncategorized")
            & (outliers["summary_comparison"] != "InsufficientData")
        ]

    # merge adjacent outliers
    outliers = merge_adjacent_outliers(outliers)
    # add a maximum delta column
    outliers["max_abs_meth_delta_zscore"] = outliers.apply(
        lambda x: (
            np.abs(x["mean_abs_meth_delta_zscore"])
            if pd.isna(x["mean_combined_methyl_zscore"])
            else (
                np.abs(x["mean_combined_methyl_zscore"])
                if pd.isna(x["mean_abs_meth_delta_zscore"])
                else np.maximum(
                    np.abs(x["mean_abs_meth_delta_zscore"]),
                    np.abs(x["mean_combined_methyl_zscore"]),
                )
            )
        ),
        axis=1,
    )
    # add a maximum number of CpGs column
    outliers["max_num_cpgs"] = outliers[
        ["num_phased_cpgs", "num_unphased_cpgs", "num_partial_cpgs"]
    ].max(axis=1)

    # convert to pyranges object
    outliers_pr = pr.PyRanges(outliers[["Chromosome", "Start", "End"]])

    # annotate with sample variants: SVs, CNVs, TRs, and small variants
    SVs_pr = SVs_to_pr(SV_path, sample)
    CNVs_pr = CNVs_to_pr(CNV_path, sample)
    TRs_pr = TRs_to_pr(TR_outlier_path, sample)
    SMVs_pr = small_variants_to_pr(SMV_path)

    outliers_sv = find_closest_variant(outliers_pr, SVs_pr, sample, "SV")
    outliers_cnv = find_closest_variant(outliers_pr, CNVs_pr, sample, "CNV")
    outliers_tr = find_closest_variant(outliers_pr, TRs_pr, sample, "TR")
    outliers_smv = find_closest_variant(outliers_pr, SMVs_pr, sample, "SMV")

    outliers_variants = (
        outliers.merge(outliers_sv, on=["Chromosome", "Start", "End"], how="left")
        .merge(outliers_cnv, on=["Chromosome", "Start", "End"], how="left")
        .merge(outliers_tr, on=["Chromosome", "Start", "End"], how="left")
        .merge(outliers_smv, on=["Chromosome", "Start", "End"], how="left")
    )

    outliers_variants["nearby_variant"] = outliers_variants.apply(
        lambda x: (";").join(
            [
                v
                for v in [
                    x["closest_SV"],
                    x["closest_CNV"],
                    x["closest_TR"],
                    x["closest_SMV"],
                ]
                if not pd.isna(v)
            ]
        ),
        axis=1,
    )

    # annotate outliers with Ensembl genes
    print("Annotate against Ensembl genes")
    outliers_pr = pr.PyRanges(outliers_variants)
    gene_gr = pd.read_csv(ensembl)
    outliers_gene = annotate_genes(outliers_pr, pr.PyRanges(gene_gr))

    # annotate with closest exon
    print("Annotate with closest exon")
    exons = gene_gr[gene_gr["Feature"] == "exon"]
    exons_pr = pr.PyRanges(exons)
    outliers_exon = find_closest_exon(outliers_gene, exons_pr)
    outliers_gene = outliers_gene.merge(
        outliers_exon, on=["Chromosome", "Start", "End"], how="left"
    )

    # annotate with gene constraint
    print("Add gnomAD gene constraint")
    constraint_cols = ["gene", "mane_select", "lof.oe_ci.upper", "lof.pLI"]
    constraint = pd.read_csv(constraint, sep="\t")[constraint_cols].dropna()
    outliers_gene = add_constraint(constraint, outliers_gene)

    # annotate with OMIM
    print("Add OMIM phenotype")
    omim = prepare_OMIM(f"{omim}/genemap2.txt")
    outliers_gene_omim = annotate_OMIM(outliers_gene, omim)

    # annotate with HPO terms
    try:
        print("Add HPO terms")
        hpo = pd.read_csv(hpo, sep="\t")
        outliers_gene_omim = add_hpo(hpo, outliers_gene_omim)
    except:
        print("No HPO terms supplied")
        outliers_gene_omim["HPO"] = ""

    # add a trid column so dataframe is compatible with group_by_gene and group_by_greendb functions
    outliers_gene_omim["trid"] = outliers_gene_omim.apply(
        lambda x: x["Chromosome"] + str(x["Start"]) + str(x["End"]), axis=1
    )
    outliers_gene_omim = group_by_gene(outliers_gene_omim)
    print(outliers_gene_omim.columns)

    # annotate with GREENDB
    print("Add GREENDB regulatory regions")
    outliers_gene_omim = annotate_reg_regions(
        outliers_gene_omim,
        "/hpf/largeprojects/ccmbio/nhanafi/c4r/downloads/databases/GRCh38_GREEN-DB.bed.gz",
    )
    outliers_gene_omim = group_by_greendb(outliers_gene_omim)

    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        outliers_gene_omim[col] = outliers_gene_omim[col].apply(
            lambda genes: gene_set(genes)
        )

    outliers_gene_omim = outliers_gene_omim.rename(
        columns={
            "Feature": "feature",
            "Segdup": "segdup",
            "Chromosome": "CHROM",
            "Start": "POS",
            "End": "END",
            "gene": "gnomad_constraint_gene",
            "Count": "tile_count",
        }
    )

    # round numeric columns
    outliers_gene_omim["asm_fishers_pvalue"] = outliers_gene_omim[
        "asm_fishers_pvalue"
    ].round(7)
    # these numeric columns have some missing values encoded as "."
    for col in [
        "mean_hap1_methyl",
        "mean_hap2_methyl",
        "mean_meth_delta",
        "mean_abs_meth_delta_zscore",
    ]:
        outliers_gene_omim[col] = outliers_gene_omim[col].replace(".", np.nan)
    numeric_cols = [
        "mean_coverage",
        "category_pop_freq",
        "mean_hap1_methyl",
        "mean_hap2_methyl",
        "mean_meth_delta",
        "mean_abs_meth_delta_zscore",
        "mean_combined_methyl",
        "mean_combined_methyl_zscore",
        "max_abs_meth_delta_zscore",
    ]
    outliers_gene_omim[numeric_cols] = outliers_gene_omim[numeric_cols].round(2)

    # Replace various missing value indicators with periods
    outliers_gene_omim.replace(
        {
            "-1": ".",
            "nan": ".",
            np.nan: ".",
            "-1;-1": ".",
            "-1;": "",
            ";-1": "",
            "": ".",
        },
        inplace=True,
    )

    columns = [
        "CHROM",
        "POS",
        "END",
        "tile_count",
        "summary_label",
        "compare_label",
        "gene_name",
        "gene_id",
        "gene_biotype",
        "closest_exon_gene",
        "omim_phenotype",
        "omim_inheritance",
        "HPO",
        "gnomad_constraint_gene",
        "lof.oe_ci.upper",
        "lof.pLI",
        "feature",
        "nearby_variant",
        "category_pop_count",
        "category_pop_freq",
        "asm_fishers_pvalue",
        "mean_hap1_methyl",
        "mean_hap2_methyl",
        "mean_meth_delta",
        "mean_abs_meth_delta_zscore",
        "mean_combined_methyl",
        "mean_combined_methyl_zscore",
        "max_abs_meth_delta_zscore",
        "num_phased_cpgs",
        "num_partial_cpgs",
        "num_unphased_cpgs",
        "max_num_cpgs",
        "mean_coverage",
        "GREENDB_reg_region",
        "GREENDB_source",
        "GREENDB_closest_gene",
        "GREENDB_controlled_genes",
    ]

    try:
        outliers_gene_omim = outliers_gene_omim[columns]
        outliers_gene_omim = outliers_gene_omim.sort_values(
            by="max_abs_meth_delta_zscore", ascending=False
        )
    except KeyError:
        columns = [
            "CHROM",
            "POS",
            "END",
            "tile_count",
            "summary_label",
            "compare_label",
            "gene_name",
            "gene_id",
            "gene_biotype",
            "omim_phenotype",
            "omim_inheritance",
            "gnomad_constraint_gene",
            "lof.oe_ci.upper",
            "lof.pLI",
            "feature",
            "nearby_variant",
            "category_pop_count",
            "category_pop_freq",
            "asm_fishers_pvalue",
            "mean_hap1_methyl",
            "mean_hap2_methyl",
            "mean_meth_delta",
            "mean_abs_meth_delta_zscore",
            "mean_combined_methyl",
            "mean_combined_methyl_zscore",
            "max_abs_meth_delta_zscore",
            "num_phased_cpgs",
            "num_partial_cpgs",
            "num_unphased_cpgs",
            "max_num_cpgs",
            "GREENDB_reg_region",
            "GREENDB_source",
            "GREENDB_closest_gene",
            "GREENDB_controlled_genes",
        ]
        outliers_gene_omim = outliers_gene_omim[columns]
        outliers_gene_omim = outliers_gene_omim.sort_values(
            by="max_abs_meth_delta_zscore", ascending=False
        )

    # write to file
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    out_file = out_file.replace(".csv", "")
    outliers_gene_omim = outliers_gene_omim.drop_duplicates()
    print(f"Writing to file {out_file}.{today}.csv")
    outliers_gene_omim.to_csv(f"{out_file}.csv", index=False)
    outliers_gene_omim.to_csv(f"{out_file}.{today}.csv", index=False)


if __name__ == "__main__":
    # if running from the command-line
    description = "Annotate repeat outliers with Ensembl genes, gnomAD constraint metrics, OMIM phenotypes, and HPO terms"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--outliers",
        type=str,
        required=True,
        help="Methylation outliers called by MethBat",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )
    parser.add_argument(
        "--ensembl",
        type=str,
        required=True,
        help="Path to Ensembl gene CSV",
    )
    parser.add_argument(
        "--gnomad_constraint",
        type=str,
        required=True,
        help="Path to gnomAD constraint file",
    )
    parser.add_argument(
        "--OMIM_path",
        type=str,
        required=True,
        help="Path to directory containing OMIM mim2gene and morbidmap files",
    )
    parser.add_argument(
        "--SV_path",
        type=str,
        required=True,
        help="Path to SV report",
    )
    parser.add_argument(
        "--CNV_path",
        type=str,
        required=True,
        help="Path to TCAG CNV report",
    )
    parser.add_argument(
        "--TR_outlier_path",
        type=str,
        required=True,
        help="Path to TR outlier report",
    )
    parser.add_argument(
        "--SMV_path",
        type=str,
        required=True,
        help="Path to annotated small variant VCF",
    )
    parser.add_argument(
        "--coverage",
        type=str,
        help="mosdepth regions.bed.gz file for 200bp tiles",
    )
    parser.add_argument(
        "--hpo",
        type=str,
        help="Path to HPO terms file",
    )

    args = parser.parse_args()
    print("Annotating methylation outliers...")
    if args.coverage:
        main(
            args.outliers,
            args.output_file,
            args.ensembl,
            args.gnomad_constraint,
            args.OMIM_path,
            args.SV_path,
            args.CNV_path,
            args.TR_outlier_path,
            args.SMV_path,
            args.coverage,
            args.hpo,
        )
    else:
        print("Cohort comparison provided, omitting HPO and coverage annotations")
        main(
            args.outliers,
            args.output_file,
            args.ensembl,
            args.gnomad_constraint,
            args.OMIM_path,
            args.SV_path,
            args.CNV_path,
            args.TR_outlier_path,
            args.SMV_path,
        )
