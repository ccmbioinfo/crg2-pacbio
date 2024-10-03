import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr
from typing import Optional

from annotation.annotate import (
    load_sample,
    hits_to_pr,
    filter_candidates,
    annotate_genes,
    gene_set,
    add_constraint,
    prepare_OMIM,
    annotate_OMIM,
    add_hpo,
    group_by_gene, 
    compare_to_controls,
    get_denovo_al,
    annotate_segdup,
    group_by_segdup
)



def main(
    hits: str,
    out_file: str,
    ensembl: str,
    constraint: str,
    omim: str,
    controls: str,
    segdup: str,
    hpo: Optional[str] = None,
    c4r: Optional[str] = None
) -> None:

    print("Loading and filtering TRGT-denovo output")
    hits = load_sample(hits)
    # filter on TRGT-denovo metrics 
    hits = filter_candidates(hits)
    # convert hits to PyRanges object
    hits_pr = hits_to_pr(hits)

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

    # annotate with segmental duplications
    segdup = pd.read_csv(segdup, sep="\t", compression="gzip", header=None, names=["Chromosome", "Start", "End", "Segdup", "Strand", "Score"])
    segdup["Chromosome"] = segdup["Chromosome"].str.replace("chr", "")
    hits_gene_omim = annotate_segdup(pr.PyRanges(hits_gene_omim), pr.PyRanges(segdup))
    hits_gene_omim = group_by_segdup(hits_gene_omim)

    # annotate with control allele distributions
    print("Annotate against control distributions")
    controls = pd.read_csv(controls, sep="\t")
    controls.drop(columns=["MP_mean", "MP_std", "AM_mean", "AM_std"], inplace=True)
    # controls["merge_key"] = controls["trid"].str.rsplit("_", n=1).str[0] # CMH TRIDs don't always have motif? Possible issue with mixed TRGT versions
    # controls.drop(columns=["trid"], inplace=True)
    controls.rename(columns={"trid": "merge_key"}, inplace=True)
    hits_gene_omim["merge_key"] = hits_gene_omim["trid"].str.rsplit("_", n=1).str[0]
    hits_gene_omim = hits_gene_omim.merge(controls, on="merge_key", how="left")

    # split into short and long alleles
    print("Split into short and long alleles") 
    short_alleles = hits_gene_omim[hits_gene_omim["genotype"] == 1].copy()
    short_alleles.drop(columns=[col for col in short_alleles.columns if "long" in col], inplace=True)
    short_alleles.rename(columns={"cutoff_short": "cutoff", "range_short": "control_range", 
                                  "short_allele_len_mean": "allele_len_mean", 
                                  "short_allele_len_std": "allele_len_std"}, inplace=True)

    long_alleles = hits_gene_omim[hits_gene_omim["genotype"] == 2].copy()
    long_alleles.drop(columns=[col for col in long_alleles.columns if "short" in col], inplace=True)
    long_alleles.rename(columns={"cutoff_long": "cutoff", "range_long": "control_range", 
                                 "long_allele_len_mean": "allele_len_mean", 
                                 "long_allele_len_std": "allele_len_std"}, inplace=True)

    # concatenate short and long alleles
    hits_gene_omim = pd.concat([short_alleles, long_alleles])

    # get length of de novo allele
    hits_gene_omim["child_denovo_AL"] = hits_gene_omim.apply(lambda x: get_denovo_al(x["child_AL"], x["genotype"]), axis=1)

    # calculate allele length z scores relative to controls and whether or not allele is an outlier
    print("Determining outliers") 
    hits_gene_omim['z_score_len'], hits_gene_omim['outlier'] = zip(*hits_gene_omim.apply(lambda x: compare_to_controls(x['child_AL'], x['genotype'], x['allele_len_mean'], x['allele_len_std'], x['cutoff']), axis=1))

    # add allele type: 1 is short allele, 2 is long allele
    hits_gene_omim["genotype"] = hits_gene_omim["genotype"].apply(lambda x: "long_allele" if x == 2 else "short_allele")
    hits_gene_omim["trid"] = hits_gene_omim["trid"].astype(str) + "_" + hits_gene_omim["genotype"].astype(str)
    hits_gene_omim.rename({"trid": "TRID"}, axis=1, inplace=True)

    # annotate with C4R outlier counts
    print("Adding inhouse outlier counts")
     
    try: 
        c4r_counts = pd.read_csv(c4r)
        hits_gene_omim = hits_gene_omim.merge(c4r_counts, on="TRID", how="left")
        hits_gene_omim["count"] = hits_gene_omim["count"].replace(np.nan, 0)
        hits_gene_omim = hits_gene_omim.rename({"count": "C4R_outlier_count", "samples": "C4R_outlier_samples"}, axis=1)
        c4r_col = ["C4R_outlier_count", "C4R_outlier_samples"]
    except:
        print("Failed to add C4R counts") 
        c4r_col = []


    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        hits_gene_omim[col] = hits_gene_omim[col].apply(lambda genes: gene_set(genes))

    # after filters, all dropout columns are 'NA', *overlap_coverage is 0
    # we can drop the genotype columns as it indicates only if the alleles is short or long, which is now integrated into the TRID
    drop_cols = [col for col in hits_gene_omim.columns if 'dropout' in col or 'index' in col or 'overlap_coverage' in col or 'genotype' in col]
    hits_gene_omim.drop(columns=drop_cols, inplace=True, errors='ignore')
    # control distribution columns
    control_distribution_columns = [
        "control_range",
        "cutoff",
    ]

    # order columns
    trgt_denovo_al_cols = ["child_AL", "child_denovo_AL", "mother_AL", "father_AL"]
    trgt_denovo_cols = [
        "denovo_coverage", "allele_coverage", "allele_ratio", "child_coverage", 
        "child_ratio", "mean_diff_father", "mean_diff_mother", "allele_origin", 
        "denovo_status", "per_allele_reads_father", "per_allele_reads_mother", 
        "per_allele_reads_child", "father_MC", "mother_MC", "child_MC", "min_mean_diff"
    ]
    hits_gene_omim = hits_gene_omim[
        ["Chromosome", "Start", "End", "TRID", "gene_name", "gene_id", "gene_biotype", "Segdup"]
        + ["omim_phenotype", "omim_inheritance", "HPO"]
        + constraint_cols
        + ["Feature"]
        + trgt_denovo_al_cols
        + ["outlier"]
        + ["z_score_len"]
        + control_distribution_columns
        + c4r_col
        + trgt_denovo_cols
    ]

    hits_gene_omim = hits_gene_omim.rename(
        columns={
            "Feature": "feature",
            "Segdup": "segdup",
            "Chromosome": "CHROM",
            "Start": "POS",
            "End": "END",
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
    hits_gene_omim.replace({"-1": ".", "nan": "."}, inplace=True)

    # write to file
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    out_file = out_file.replace(".csv", "")
    hits_gene_omim = hits_gene_omim.drop_duplicates()
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
        "--controls",
        type=str,
        required=True,
        help="Control allele length distribution TSV",
    )
    parser.add_argument(
        "--segdup",
        type=str,
        required=True,
        help="Segmental duplication BED",
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
        args.ensembl,
        args.gnomad_constraint,
        args.OMIM_path,
        args.controls,
        args.segdup,
        args.hpo,
        args.c4r_outliers
    )
