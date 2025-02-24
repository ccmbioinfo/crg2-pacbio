import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr
import sys
from typing import Optional

sys.path.insert(0, "/home/ccmmarvin/crg2-pacbio/scripts/annotation")

def annotate_adjacents_tiles(chr: str, start: int, end: int, coordinates_dict: dict) -> int:
    """determine if outlier region is adjacent to 0, one, or two outlier regions"""
    count = 0
    start = f"{chr}-{start}"
    end = f"{chr}-{end}"
    if coordinates_dict[start] > 1: 
        count += 1
    if coordinates_dict[end] > 1:
        count += 1

    return count 

def find_closest_exon(outlier_row: pd.Series, exons: pd.DataFrame) -> pd.Series:
    """Find the closest exon to an outlier region"""
    # Filter exons on same chromosome
    chrom_exons = exons[exons["Chromosome"] == outlier_row["Chromosome"]]
    
    if len(chrom_exons) == 0:
        return pd.Series({
            "closest_exon_dist": None,
            "closest_exon_gene": None,
            "closest_exon_start": None,
            "closest_exon_end": None
        })
    
    # Vectorized distance calculation: calculate start and end distances for all exons at once
    # Where the start of the exon is greater than the end of the outlier, the distance is the difference between the start of the exon and the end of the outlier
    # Where the end of the exon is less than the start of the outlier, the distance is the difference between the end of the exon and the start of the outlier
    # Where the exon overlaps with the outlier, the distance is 0
    start_dists = np.where(outlier_row["End"] < chrom_exons["Start"],
                          chrom_exons["Start"] - outlier_row["End"],
                          np.inf)
    
    end_dists = np.where(outlier_row["Start"] > chrom_exons["End"], 
                         outlier_row["Start"] - chrom_exons["End"],
                         np.inf)
    
    # Check for overlaps between outlier region and exons
    # overlaps is a series of booleans
    overlaps = ((outlier_row["Start"] <= chrom_exons["End"]) & 
                (outlier_row["End"] >= chrom_exons["Start"]))
    
    # Combine all distances, using 0 for overlaps
    distances = np.where(overlaps, 0, 
                        np.minimum(start_dists, end_dists))
    
    # Find closest exon
    min_dist_idx = np.argmin(distances)
    closest_exon = chrom_exons.iloc[min_dist_idx]
    
    return pd.Series({
        "closest_exon_dist": distances[min_dist_idx],
        "closest_exon_gene": closest_exon["gene_name"] if pd.isna(closest_exon["gene_name"]) == False else closest_exon["gene_id"],
        "closest_exon_start": closest_exon["Start"],
        "closest_exon_end": closest_exon["End"]
    })

def annotate_reg_regions(outliers: pd.DataFrame, greendb: str) -> pd.DataFrame:
    """Annotate outliers with GREENDB regulatory regions"""
    greendb = pd.read_csv(greendb, sep="\t", compression="gzip")
    greendb.rename(columns={"#Chromosome": "Chromosome"}, inplace=True)
    # convert greendb and outlier dataframes to PyRanges objects and join them to find overlaps
    greendb_pr = pr.PyRanges(greendb)
    outliers.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}, inplace=True)
    outliers_pr = pr.PyRanges(outliers)
    outliers_pr = outliers_pr.join(greendb_pr, how="left", suffix="_greendb")
    outliers = outliers_pr.df.drop(columns=["Start_greendb", "End_greendb", "regionID", "constrain_pct", "PhyloP100_median", "closestGene_dist", "closestProt_symbol", "closestProt_dist", "N_methods"])
    outliers.rename(columns={"std_type": "GREENDB_reg_region", "db_source": "GREENDB_source", "closestGene_symbol": "GREENDB_closest_gene", "controlled_genes": "GREENDB_controlled_genes"}, inplace=True)
    outliers.replace({np.nan: "."}, inplace=True)

    return outliers

def group_by_greendb(hits_gene: pd.DataFrame) -> pd.DataFrame:
    """
    One hit may be associated with multiple GREENDB features and therefore multiple rows
    Aggregate by GREENDB and join features to remove duplicate rows
    """
    hits_greendb_dedup = hits_gene.groupby(["trid"]).agg(
        {
            "GREENDB_reg_region": ";".join,
            "GREENDB_source": ";".join,
            "GREENDB_closest_gene": ";".join,
            "GREENDB_controlled_genes": ";".join,
        }
    )
    # merge with original loci table
    hits_gene = hits_gene.drop(columns=["GREENDB_reg_region", "GREENDB_source", "GREENDB_closest_gene", "GREENDB_controlled_genes"])
    hits_greendb_merged = hits_gene.merge(hits_greendb_dedup, on=["trid"], how="left")
    hits_greendb_merged_dedup = hits_greendb_merged.drop_duplicates(keep="first")

    return hits_greendb_merged_dedup
from annotate import (
    annotate_genes,
    gene_set,
    add_constraint,
    prepare_OMIM,
    annotate_OMIM,
    add_hpo,
    group_by_gene,
)

from annotate import (
    annotate_genes,
    gene_set,
    add_constraint,
    prepare_OMIM,
    annotate_OMIM,
    add_hpo,
    group_by_gene,
)

def main(
    outliers: str,
    out_file: str,
    ensembl: str,
    constraint: str,
    omim: str,
    coverage: Optional[str] = None,
    hpo: Optional[str] = None
) -> None:

    print("Loading and filtering methylation outliers")
    outliers = pd.read_csv(outliers, sep="\t")
    try:
        coverage = pd.read_csv(coverage, compression="gzip", sep="\t", header=None, names=["chrom", "start", "end", "cpg_ID", "mean_coverage"])
        # annotate regions with average coverage from mosdepth
        outliers = outliers.merge(coverage, on=["chrom", "start", "end"], how="left")
    except:
        print("No coverage file provided")
  
    try:
        outliers["chrom"] = outliers["chrom"].str.replace("chr", "")
    except KeyError: #cohort comparison
        outliers.rename(columns={"#chrom": "chrom"}, inplace=True)
        outliers["chrom"] = outliers["chrom"].str.replace("chr", "")

    try:
        outliers = outliers[(outliers["compare_label"] != "Uncategorized") & 
                            (outliers["compare_label"] != "InsufficientData")]
    except KeyError: # cohort comparison
        outliers = outliers[(outliers["summary_comparison"] != "Uncategorized") & 
                            (outliers["summary_comparison"] != "InsufficientData")]

    # count number of adjacent outlier tiles
    outliers["chr-start"] = outliers["chrom"] + "-" + outliers["start"].astype(str)
    outliers["chr-end"] = outliers["chrom"] + "-" + outliers["end"].astype(str)
    coordinates = outliers["chr-start"].tolist()
    coordinates.extend(outliers["chr-end"].tolist())
    coordinates_dict = {}
    for coord in coordinates:
        coordinates_dict[coord] = coordinates.count(coord)
    outliers["adjacent_tiles"] = outliers.apply(lambda row: annotate_adjacents_tiles(row["chrom"], row["start"], row["end"], coordinates_dict), axis=1)

    # convert to pyranges object
    outliers.rename({"chrom": "Chromosome", "start": "Start", "end": "End"}, axis=1, inplace=True)
    outliers_pr = pr.PyRanges(outliers)

    # annotate outliers with Ensembl genes
    print("Annotate against Ensembl genes")
    gene_gr = pd.read_csv(ensembl)
    outliers_gene = annotate_genes(outliers_pr, pr.PyRanges(gene_gr))

    # annotate with closest exon
    print("Annotate with closest exon")
    exons = gene_gr[gene_gr["Feature"] == "exon"]
    closest_exons = outliers_gene.apply(lambda x: find_closest_exon(x, exons), axis=1)
    outliers_gene = pd.concat([outliers_gene, closest_exons], axis=1)


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
    outliers_gene_omim["trid"] = outliers_gene_omim.apply(lambda x: x["Chromosome"] + str(x["Start"]) + str(x["End"]), axis=1)
    print(outliers_gene_omim.columns)
    outliers_gene_omim = group_by_gene(outliers_gene_omim)

    # annotate with GREENDB
    print("Add GREENDB regulatory regions")
    outliers_gene_omim = annotate_reg_regions(outliers_gene_omim, "/hpf/largeprojects/ccmbio/nhanafi/c4r/downloads/databases/GRCh38_GREEN-DB.bed.gz")
    print(outliers_gene_omim.columns)
    outliers_gene_omim = group_by_greendb(outliers_gene_omim)


    # add a maximum delta column
    outliers_gene_omim['max_abs_meth_delta_zscore'] = outliers_gene_omim.apply(
        lambda x: np.abs(x['mean_abs_meth_delta_zscore']) if pd.isna(x['mean_combined_methyl_zscore'])
        else np.abs(x['mean_combined_methyl_zscore']) if pd.isna(x['mean_abs_meth_delta_zscore'])
        else np.maximum(np.abs(x['mean_abs_meth_delta_zscore']), np.abs(x['mean_combined_methyl_zscore'])),
        axis=1
    )

    # add a maximum number of CpGs column
    outliers_gene_omim['max_num_cpgs'] = outliers_gene_omim[['num_phased_cpgs', 'num_unphased_cpgs']].max(axis=1)
    
    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        outliers_gene_omim[col] = outliers_gene_omim[col].apply(lambda genes: gene_set(genes))
    

    outliers_gene_omim = outliers_gene_omim.rename(
        columns={
            "Feature": "feature",
            "Segdup": "segdup",
            "Chromosome": "CHROM",
            "Start": "POS",
            "End": "END",
            "gene": "gnomad_constraint_gene",
        }
    )

    outliers_gene_omim.replace({"-1": ".", "nan": "."}, inplace=True)


    columns = [
        "CHROM", "POS", "END", "adjacent_tiles", "summary_label", "compare_label", 
        "gene_name", "gene_id", "gene_biotype", "closest_exon_gene", "omim_phenotype", "omim_inheritance", 
        "HPO", "gnomad_constraint_gene", "lof.oe_ci.upper", "lof.pLI", "feature",
        "category_pop_count", "category_pop_freq", "asm_fishers_pvalue",
        "mean_hap1_methyl", "mean_hap2_methyl", "mean_meth_delta",
        "mean_abs_meth_delta_zscore", "mean_combined_methyl",
        "mean_combined_methyl_zscore", "max_abs_meth_delta_zscore", "num_phased_cpgs", "num_partial_cpgs",
        "num_unphased_cpgs", "max_num_cpgs", "mean_coverage", "GREENDB_reg_region", "GREENDB_source", "GREENDB_closest_gene", "GREENDB_controlled_genes"
    ]

    try:
        outliers_gene_omim = outliers_gene_omim[columns]
        outliers_gene_omim = outliers_gene_omim.sort_values(by="max_abs_meth_delta_zscore", ascending=False)
    except KeyError:
    #     columns = [
    #         "CHROM", "POS", "END", "baseline_category", "compare_category", "summary_comparison",
    #         "gene_name", "gene_id", "gene_biotype", "omim_phenotype", "omim_inheritance", 
    #         "HPO", "gnomad_constraint_gene", "lof.oe_ci.upper", "lof.pLI", "feature",
    #         "zscore_avg_abs_meth_deltas",
    #         "delta_avg_abs_meth_deltas",
    #         "baseline_num_phased",
    #         "compare_num_phased",
    #         "zscore_avg_combined_methyls",
    #         "delta_avg_combined_methyls", 
    #         "baseline_num_samples",
    #         "compare_num_samples"
    # ]
        columns = [
            "CHROM", "POS", "END", "summary_label", "compare_label", 
            "gene_name", "gene_id", "gene_biotype", "omim_phenotype", "omim_inheritance",  "gnomad_constraint_gene", "lof.oe_ci.upper", "lof.pLI", "feature",
            "category_pop_count", "category_pop_freq", "asm_fishers_pvalue",
            "mean_hap1_methyl", "mean_hap2_methyl", "mean_meth_delta",
            "mean_abs_meth_delta_zscore", "mean_combined_methyl",
            "mean_combined_methyl_zscore", "max_abs_meth_delta_zscore", "num_phased_cpgs", "num_partial_cpgs",
            "num_unphased_cpgs", "max_num_cpgs",  "GREENDB_reg_region", "GREENDB_source", "GREENDB_closest_gene", "GREENDB_controlled_genes"
        ]
        outliers_gene_omim = outliers_gene_omim[columns]
        outliers_gene_omim = outliers_gene_omim.sort_values(by="max_abs_meth_delta_zscore", ascending=False)


    # write to file
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    out_file = out_file.replace(".csv", "")
    outliers_gene_omim = outliers_gene_omim.drop_duplicates()
    print(f"Writing to file {out_file}.{today}.csv")
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
        "--coverage",
        type=str,
        help="mosdepth regions.bed.gz file for CpG islands",
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
        )




