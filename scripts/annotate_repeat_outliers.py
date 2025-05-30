import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr
import re
from typing import Optional

from annotation.annotate import (
    pivot_hits,
    filter_outliers,
    hits_to_pr,
    annotate_motif,
    annotate_genes,
    gene_set,
    num_expanded,
    add_constraint,
    prepare_OMIM,
    annotate_OMIM,
    add_hpo,
    group_by_gene,
    calculate_reciprocal_overlap
)

def main(
    hits: pd.DataFrame,
    out_file: str,
    ensembl: str,
    constraint: str,
    omim: str,
    promoters: str,
    TR_constraint: str,
    c4r: bool,
    repeat_catalog: str,
    hpo: Optional[str] = None,
    c4r_counts: Optional[str] = None
) -> None:
    # convert hits dataframe from long to wide format
    hits = pd.read_csv(hits, sep="\t")
    hits_pivot = pivot_hits(hits)

    # make a column with maximum z score for allele length across samples
    z_score_cols = [col for col in hits_pivot.columns if "z_score_len" in col and "rank" not in col]
    z_score_len_rank_cols = [col for col in hits_pivot.columns if "z_score_len_rank" in col]
    for col in z_score_cols:
        hits_pivot[col] = [
            round(score, 3) if not pd.isnull(score) else None
            for score in hits_pivot[col]
        ]
    # filter out non-outliers
    print("Filter outliers by outlier cutoff and Z-score")
    al_cols = [col for col in hits_pivot.columns if "_allele_len" in col]
    hits_pivot["outlier"] = hits_pivot.apply(
        lambda row: filter_outliers(row, al_cols), axis=1
    )
    hits_pivot = hits_pivot[hits_pivot["outlier"]]
    hits_pivot["max_z_score_len"] = hits_pivot[z_score_cols].max(axis=1)
    hits_pivot = hits_pivot[hits_pivot["max_z_score_len"] >= 3]


    hits_pivot["min_z_score_len_rank"] = hits_pivot[z_score_len_rank_cols].min(axis=1)

    # make a column with maximum LPS across samples
    lps_cols = [col for col in hits_pivot.columns if "LPS" in col and 'rank' not in col and 'z_score' not in col]
    for col in lps_cols:
        hits_pivot[col] = [
            lps if not pd.isnull(lps) else None
            for lps in hits_pivot[col]
        ]
    hits_pivot["max_LPS"] = hits_pivot[lps_cols].max(axis=1)

    # make a column with maximum allele length across samples
    al_cols = [col for col in hits_pivot.columns if "_allele_len" in col]
    hits_pivot["max_sample_allele_len"] = hits_pivot[al_cols].max(axis=1)

    # annotate with motif
    print("Annotate with motif")
    hits_pivot = annotate_motif(hits_pivot, repeat_catalog)

    # make a column that sums the number of individuals carrying a particular repeat expansion
    hits_pivot["num_outlier_samples"] = hits_pivot.apply(
        lambda row: num_expanded(row, al_cols, z_score_cols), axis=1
    )

    # convert hits to PyRanges object
    hits_pr = hits_to_pr(hits_pivot)

    # annotate hits with Ensembl genes
    print("Annotate against Ensembl genes")
    gene_gr = pd.read_csv(ensembl)
    hits_gene = annotate_genes(hits_pr, pr.PyRanges(gene_gr))

    # annotate with gene constraint
    print("Add gnomAD gene constraint")
    constraint_cols = ["gene", "mane_select", "lof.oe_ci.upper", "lof.pLI"]
    constraint = pd.read_csv(constraint, sep="\t")[constraint_cols].dropna()
    hits_gene = add_constraint(constraint, hits_gene)

    # annotate with OMIM
    print("Add OMIM phenotype")
    omim = prepare_OMIM(f"{omim}/genemap2.txt")
    hits_gene_omim = annotate_OMIM(hits_gene, omim)

    # annotate with HPO terms
    if hpo == None:
        print("No HPO terms supplied")
        hits_gene_omim["HPO"] = ""
    else:
        print("Add HPO terms")
        hpo = pd.read_csv(hpo, sep="\t")
        hits_gene_omim = add_hpo(hpo, hits_gene_omim)

    # group and aggregate gene columns
    hits_gene_omim = group_by_gene(hits_gene_omim)

    # annotate with C4R outlier counts
    print("Adding inhouse outlier counts")
    print(c4r)
     
    try: 
        c4r_counts = pd.read_csv(c4r_counts)
        c4r_counts["TRID"] = c4r_counts["TRID"].str.rsplit("_", n=3).str[0] + "_" + c4r_counts["TRID"].str.split("_").str[4]
        hits_gene_omim = hits_gene_omim.merge(c4r_counts, left_on="trid", right_on="TRID", how="left")
        hits_gene_omim["count"] = hits_gene_omim["count"].replace(np.nan, 0)
        hits_gene_omim = hits_gene_omim.rename({"count": "C4R_outlier_count", "samples": "C4R_outlier_samples", "allele_lens": "C4R_outlier_allele_lens"}, axis=1)
        c4r_col = ["C4R_outlier_count", "C4R_outlier_samples", "C4R_outlier_allele_lens"]
    except:
        print("Failed to add C4R counts") 
        c4r_col = []
    
    if c4r != "True": 
        c4r_col = ["C4R_outlier_count", "C4R_outlier_allele_lens"] # remove C4R sample IDs for non-C4R projects


    # add promoters
    promoters = pr.read_bed(promoters)
    hits_gene_omim_pr = pr.PyRanges(hits_gene_omim)
    hits_gene_omim = hits_gene_omim_pr.join(promoters, how="left", suffix="_promoter").df
    hits_gene_omim.loc[hits_gene_omim["Score"] == "-1", "ENCODE_promoter_coord"] = "."
    hits_gene_omim.loc[hits_gene_omim["Score"] != "-1", "ENCODE_promoter_coord"] = hits_gene_omim["Chromosome"].astype(str) + ":" + hits_gene_omim["Start_promoter"].astype(str) + "-" + hits_gene_omim["End_promoter"].astype(str)
    hits_gene_omim.rename(columns={"Score": "ENCODE_promoter_ID"}, inplace=True)
    hits_gene_omim = hits_gene_omim.drop(columns=["Start_promoter", "End_promoter", "Name", "Strand"])

    # add TR constraints
    # supplementary table 3 from 'Detailed tandem repeat allele profiling in 1,027 long-read genomes reveals genome-wide patterns of pathogenicity'
    constraints = pd.read_csv(TR_constraint, dtype={"TRID": str, "combinedLPSStdev": float, "expectedCombinedLPSStdev": float, "OE_len": float, "CPS": float, "expectedCPS": float, "OE_motif": float})
    constraints["OE_len"] = constraints["OE_len"].round(3)
    hits_gene_omim["TRID_trim"] = hits_gene_omim["TRID"].str.rsplit("_", n=3).str[0]
    constraint_bed = pd.DataFrame.from_dict({'Chromosome': constraints["TRID"].str.split("_").str[0], 
                                            'Start': constraints["TRID"].str.split("_").str[1], 
                                            'End': constraints["TRID"].str.split("_").str[2], 
                                            'OE_len': constraints["OE_len"]})
    
    constraint_bed["Chromosome"] = constraint_bed["Chromosome"].str.replace("chr", "")
    constraint_pr = pr.PyRanges(constraint_bed)
    hits_gene_omim_pr = pr.PyRanges(hits_gene_omim)
    hits_gene_omim = hits_gene_omim_pr.join(constraint_pr, how="left", suffix="_constraint", report_overlap=True).df
    hits_gene_omim.replace({-1: np.nan}, inplace=True)

    hits_gene_omim["len_constraint"] = hits_gene_omim["End_constraint"] - hits_gene_omim["Start_constraint"]
    hits_gene_omim["len_TRID"] = hits_gene_omim["End"] - hits_gene_omim["Start"]
    hits_gene_omim["len_ratio"] = hits_gene_omim["len_constraint"] / hits_gene_omim["len_TRID"]

    # calculate reciprocal overlap of TR locus and constraint region
    hits_gene_omim['reciprocal_overlap'] = hits_gene_omim.apply(calculate_reciprocal_overlap, axis=1)

    # replace OE_len with np.nan where reciprocal_overlap threshold (0.8) not met
    hits_gene_omim.loc[hits_gene_omim['reciprocal_overlap'] < 0.8, 'OE_len'] = np.nan
    hits_gene_omim.loc[hits_gene_omim['Overlap'] < 0, 'Overlap'] = np.nan
    hits_gene_omim.drop(columns=["Start_constraint", "End_constraint", "len_constraint", "len_TRID", "len_ratio"], inplace=True)

    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        hits_gene_omim[col] = hits_gene_omim[col].apply(lambda genes: gene_set(genes))

    am_cols = [col for col in hits_gene_omim.columns if "AM" in col]
    mp_cols = [col for col in hits_gene_omim.columns if "MP" in col]
    z_score_ranks_cols = [col for col in hits_gene_omim.columns if "z_score_len_rank" in col and 'min' not in col]
    lps_cols = [col for col in hits_gene_omim.columns if "LPS" in col and 'max' not in col]

    hits_gene_omim = hits_gene_omim[
        ["Chromosome", "Start", "End", "trid", "motif", "gene_name", "gene_id", "gene_biotype", "Feature"]
        + ["omim_phenotype","omim_inheritance", "HPO"]
        + ["gene", "lof.oe_ci.upper", "lof.pLI", "OE_len"]
        + ["ENCODE_promoter_ID", "ENCODE_promoter_coord", "range", "cutoff", "allele_len_std"]
        + c4r_col
        + ["max_sample_allele_len", "max_z_score_len", "min_z_score_len_rank", "max_LPS"]
        + lps_cols
        + ["num_outlier_samples"]
        + al_cols
        + z_score_cols
        + z_score_ranks_cols
        + am_cols
        + mp_cols
    ]
    hits_gene_omim = hits_gene_omim.rename(
        columns={
            "Feature": "feature",
            "Chromosome": "CHROM",
            "Start": "POS",
            "End": "END",
            "trid": "TRID",
            "gene": "gnomad_constraint_gene", 
            "allele_len_std": "CMH_allele_len_std",
            "range": "CMH_allele_len_range",
            "cutoff": "CMH_outlier_cutoff",
            "OE_len": "TR_constraint"
        }
    )

    # sort by max_LPS
    hits_gene_omim = hits_gene_omim.sort_values(by="max_LPS", ascending=False)

    # recode dtypes and missing values
    hits_gene_omim["CHROM"] = hits_gene_omim["CHROM"].astype(str)
    try:
    	hits_gene_omim["C4R_outlier_count"] = hits_gene_omim["C4R_outlier_count"].astype(int)
    except KeyError:
        pass 
    hits_gene_omim.fillna(".", inplace=True)
    hits_gene_omim.replace({"-1": ".", "1:-1--1": ".", " ": ".", "nan": "."}, inplace=True)

    # drop dups
    hits_gene_omim = hits_gene_omim.drop_duplicates()

    # write to file
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    out_file = out_file.replace(".csv", "")
    hits_gene_omim.to_csv(f"{out_file}.csv", index=False)
    hits_gene_omim.to_csv(f"{out_file}.{today}.csv", index=False)


if __name__ == "__main__":
    # if running from the command-line
    description = "Annotate repeat outliers with Ensembl genes, gnomAD constraint metrics, OMIM phenotypes, and HPO terms"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--repeats",
        type=str,
        required=True,
        help="Repeat outlilers from find_repeat_outliers",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )
    parser.add_argument(
        "--ensembl_gtf",
        type=str,
        required=True,
        help="Path to Ensembl gene GTF",
    )
    parser.add_argument(
        "--gnomad_constraint",
        type=str,
        required=True,
        help="Path to Ensembl gene GTF",
    )
    parser.add_argument(
        "--OMIM_path",
        type=str,
        required=True,
        help="Path to directory containing OMIM mim2gene and morbidmap files",
    )
    parser.add_argument(
        "--promoters",
        type=str,
        required=True,
        help="Path to directory containing OMIM mim2gene and morbidmap files",
    )
    parser.add_argument(
        "--TR_constraint",
        type=str,
        required=True,
        help="Path to directory containing OMIM mim2gene and morbidmap files",
    )
    parser.add_argument(
        "--c4r",
        type=str,
        required=True,
        help="True if C4R project, otherwise false (C4R outlier sample IDs will be excluded from report)",
    )
    parser.add_argument(
        "--repeat_catalog",
        type=str,
        help="Path to TRGT repeat catalog",
    )
    parser.add_argument(
        "--hpo",
        type=str,
        help="Path to HPO terms file",
    )
    parser.add_argument(
        "--c4r_outliers",
        type=str,
        help="Path to C4R tandem repeat outlier count",
    )

    args = parser.parse_args()
    print("Annotating repeats...")
    main(
        args.repeats,
        args.output_file,
        args.ensembl_gtf,
        args.gnomad_constraint,
        args.OMIM_path,
        args.promoters,
        args.TR_constraint,
        args.c4r,
        args.repeat_catalog,
        args.hpo,
        args.c4r_outliers
    )
