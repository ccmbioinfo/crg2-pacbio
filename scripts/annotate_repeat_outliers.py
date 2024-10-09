import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr
import re
from typing import Optional

from annotation.annotate import (
    pivot_hits,
    hits_to_pr,
    prepare_Ensembl_GTF,
    annotate_genes,
    gene_set,
    filter_outliers,
    num_expanded,
    add_constraint,
    prepare_OMIM,
    annotate_OMIM,
    add_hpo,
    group_by_gene
)

def main(
    hits: pd.DataFrame,
    out_file: str,
    ensembl: str,
    constraint: str,
    omim: str,
    hpo: Optional[str] = None,
    c4r: Optional[str] = None
) -> None:
    # convert hits dataframe from long to wide format
    hits = pd.read_csv(hits)
    hits.rename({"case_trid": "trid"}, axis=1, inplace=True)
    hits_pivot = pivot_hits(hits)

    # make a column with maximum z score for allele length across samples
    z_score_cols = [col for col in hits_pivot.columns if "z_score_len" in col]
    for col in z_score_cols:
        hits_pivot[col] = [
            round(score, 3) if not pd.isnull(score) else None
            for score in hits_pivot[col]
        ]
    hits_pivot["max_z_score_len"] = hits_pivot[z_score_cols].max(axis=1)

    # filter out non-outliers
    print("Filter outliers")
    al_cols = [col for col in hits_pivot.columns if "_allele_len" in col]
    hits_pivot["outlier"] = hits_pivot.apply(
        lambda row: filter_outliers(row, al_cols), axis=1
    )
    hits_pivot = hits_pivot[hits_pivot["outlier"]]

    # make a column that sums the number of individuals carrying a particular repeat expansion
    hits_pivot["num_samples"] = hits_pivot.apply(
        lambda row: num_expanded(row, al_cols), axis=1
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
     
    try: 
        c4r_counts = pd.read_csv(c4r)
        hits_gene_omim = hits_gene_omim.merge(c4r_counts, left_on="trid", right_on="TRID", how="left")
        hits_gene_omim["count"] = hits_gene_omim["count"].replace(np.nan, 0)
        hits_gene_omim = hits_gene_omim.rename({"count": "C4R_outlier_count", "samples": "C4R_outlier_samples"}, axis=1)
        c4r_col = ["C4R_outlier_count", "C4R_outlier_samples"]
    except:
        print("Failed to add C4R counts") 
        c4r_col = []


    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        hits_gene_omim[col] = hits_gene_omim[col].apply(lambda genes: gene_set(genes))

    am_cols = [col for col in hits_gene_omim.columns if "AM" in col]
    mp_cols = [col for col in hits_gene_omim.columns if "MP" in col]
    

    hits_gene_omim = hits_gene_omim[
        ["Chromosome", "Start", "End", "trid", "gene_name", "gene_id", "gene_biotype"]
        + ["omim_phenotype","omim_inheritance", "HPO"]
        + ["gene", "lof.oe_ci.upper", "lof.pLI"]
        + ["Feature", "control_range", "cutoff"]
        + c4r_col
        + ["max_z_score_len", "num_samples"]
        + al_cols
        + z_score_cols
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
        }
    )

    # recode dtypes and missing values
    hits_gene_omim["CHROM"] = hits_gene_omim["CHROM"].astype(str)
    try:
    	hits_gene_omim["C4R_outlier_count"] = hits_gene_omim["C4R_outlier_count"].astype(int)
    except KeyError:
        pass 
    hits_gene_omim.fillna(".", inplace=True)
    hits_gene_omim.replace({"-1": "."}, inplace=True)

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
        args.hpo,
        args.c4r_outliers
    )
