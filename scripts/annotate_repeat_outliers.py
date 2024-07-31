import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr
from typing import Optional

def pivot_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert hits dataframe from long to wide format so each sample is associated with two columns:
    [sample]_allele_len, and [sample]_z_score
    """
    hits = df.copy()
    hits["trid"] = hits["trid"] + "_" + hits["allele_type"]
    hits = hits.drop("allele_type", axis=1)
    hit_pivot = hits.pivot(
        index=["trid", "control_range", "cutoff"],
        columns="sample",
        values=["allele_len", "z_score"],
    ).reset_index()

    hit_pivot.columns = (
        hit_pivot.columns.get_level_values(1)
        + "_"
        + hit_pivot.columns.get_level_values(0)
    )
    hit_pivot = hit_pivot.rename(
        columns={"_trid": "trid", "_control_range": "control_range", "_cutoff": "cutoff"}
    )

    return hit_pivot


def hits_to_pr(hits: pd.DataFrame) -> pr.PyRanges:
    """
    Convert hits dataframe to PyRanges object
    """
    hits["Chromosome"] = hits["trid"].map(lambda x: x.split("_")[0].replace("chr", ""))
    hits["Start"] = hits["trid"].map(lambda x: x.split("_")[1])
    hits["End"] = hits["trid"].map(lambda x: x.split("_")[2])
    hit_pr = pr.PyRanges(hits)

    return hit_pr


def prepare_Ensembl_GTF(gtf_path: str, cols: list) -> pr.PyRanges:
    """
    Convert GTF file into a PyRanges object and clean up columns
    """
    ensembl_gr = pr.read_gtf(gtf_path)[cols]
    introns = ensembl_gr.features.introns(by="gene")
    ensembl_gr = pd.concat([ensembl_gr.df, introns.df], ignore_index=True)
    features = ["exon", "three_prime_utr", "five_prime_utr", "intron"]
    ensembl_gr = ensembl_gr[ensembl_gr.Feature.isin(features)]

    return ensembl_gr


def prepare_TRGT_loci(loci_path: str) -> pr.PyRanges:
    """
    Convert repeat loci derived from TRGTdb into a PyRanges object
    """
    loci = pd.read_csv(loci_path).rename(
        {"chrom": "Chromosome", "end": "End", "start": "Start"}, axis=1
    )
    # remove chr prefix from chromosomes for compatibility with gene table
    loci["Chromosome"] = loci["Chromosome"].apply(lambda chrom: chrom.strip("chr"))
    loci_gr = pr.PyRanges(loci)

    return loci_gr


def annotate_genes(loci: pr.PyRanges, genes: pr.PyRanges) -> pd.DataFrame:
    """
    Annotate loci against Ensembl genes
    """
    loci_ensembl = loci.join(
        genes, suffix="_ensembl", how="left", apply_strand_suffix=False
    ).drop(["Start_ensembl", "End_ensembl", "Strand"])
    loci_ensembl_df = loci_ensembl.df
    loci_ensembl_df["gene_name"] = np.where(
        loci_ensembl_df["gene_name"].isnull(), "-1", loci_ensembl_df["gene_name"]
    )

    return loci_ensembl_df


def add_constraint(constraint: pd.DataFrame, hits_gene: pd.DataFrame) -> pd.DataFrame:
    """Add transcript-specific gnomAD v4 LOEUF and pLI scores"""
    hits_gene = hits_gene.merge(constraint, left_on="gene_name", right_on="gene", how="left")

    return hits_gene
    

def group_by_gene(hits_gene: pd.DataFrame) -> pd.DataFrame:
    """
    One hit may be associated with multiple gene features and so multiple rows
    Aggregate by gene and join features to remove duplicate rows
    """
    gene_cols = ["gene_name", "gene_id", "gene_biotype", "Feature", "lof.oe_ci.upper", "lof.pLI"]

    hits_gene_dedup = hits_gene.groupby(["trid"]).agg(
        {
            "gene_name": ";".join,
            "gene_id": ";".join,
            "gene_biotype": ";".join,
            "Feature": ";".join,
            "lof.oe_ci.upper": min,
            "lof.pLI": max,
        }
    )
    # merge with original loci table
    hits_gene = hits_gene.drop(gene_cols, axis=1)
    hits_gene_merged = hits_gene.merge(hits_gene_dedup, on=["trid"], how="left")
    hits_gene_merged_dedup = hits_gene_merged.drop_duplicates()

    return hits_gene_merged_dedup


def prepare_OMIM(mim2gene_path: str, morbidmap_path: str) -> pd.DataFrame:
    """
    Merge OMIM mim2gene table with morbidmap table to associate phenotypes to genes.
    See omim.org/downloads for download access to these files.
    """
    mim2gene = pd.read_csv(
        mim2gene_path,
        sep="\t",
        comment="#",
        names=[
            "MIM_Number",
            "MIM_Entry_Type",
            "Entrez_Gene_ID",
            "HGNC_symbol",
            "Ensembl_gene_ID",
        ],
    ).set_index("MIM_Number")
    morbidmap = pd.read_csv(
        morbidmap_path,
        sep="\t",
        comment="#",
        names=["Phenotype", "Gene_Symbols", "MIM_Number", "Cyto_Location"],
    ).set_index("MIM_Number")

    gene2mim = morbidmap.join(mim2gene, how="left")
    gene2mim = gene2mim[gene2mim["MIM_Entry_Type"].str.contains("gene")].drop(
        [
            "Gene_Symbols",
            "Cyto_Location",
            "MIM_Entry_Type",
            "Entrez_Gene_ID",
            "HGNC_symbol",
        ],
        axis=1,
    )
    # remove genes that are not associated with any phenotype and rename columns
    gene2mim = gene2mim[gene2mim["Phenotype"].notna()].rename(
        {"Phenotype": "OMIM_phenotype", "Ensembl_gene_ID": "gene_id"}, axis=1
    )
    # group phenotypes by gene and aggregate phenotypes so that there is one row per gene:
    gene2mim = (
        gene2mim.groupby("gene_id")[["OMIM_phenotype"]]
        .agg({"OMIM_phenotype": (";").join})
        .reset_index()
    )

    return gene2mim


def annotate_OMIM(loci_ensembl: pd.DataFrame, omim: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate loci against OMIM phenotypes. Loci must already be annotated against Ensembl genes.
    """
    loci_ensembl_omim = loci_ensembl.merge(omim, how="left", on="gene_id")
    loci_ensembl_omim["OMIM_phenotype"] = loci_ensembl_omim["OMIM_phenotype"].fillna(
        "-1"
    )

    return loci_ensembl_omim


def gene_set(genes: str) -> str:
    """
    Reduce semi-colon separated gene identifiers to unique gene identifiers
    """
    genes = ";".join(list(set(genes.split(";"))))
    return genes


def add_hpo(hpo: pd.DataFrame, loci_ensembl: pd.DataFrame) -> list:
    """
    Add gene-based HPO terms
    """
    hpo = hpo[["Gene ID", "Features"]].copy()
    hpo.rename({"Gene ID": "gene_id", "Features": "HPO"}, axis=1, inplace=True)
    loci_ensembl_hpo = loci_ensembl.merge(hpo, how="left", on="gene_id")
    loci_ensembl_hpo["HPO"] = loci_ensembl_hpo["HPO"].fillna(
        "-1"
    )

    return loci_ensembl_hpo


def filter_outliers(row: pd.Series, allele_len_cols: list) -> bool:
    """
    Filter out non-outliers (not expanded in any individual)
    """
    allele_lens = [float(row[allele_len]) if not pd.isna(row[allele_len])  else 0 for allele_len in allele_len_cols]
    cutoff = row["cutoff"]
    if max(allele_lens) < cutoff:
        return False
    else:
        return True
    

def num_expanded(row: pd.Series, allele_len_cols: list) -> bool:
    """
    Sum number of individuals in whom repeat is expanded (i.e. allele length > cutoff)
    """
    allele_lens = [row[allele_len] for allele_len in allele_len_cols]
    cutoff = row["cutoff"] if not pd.isna(row["cutoff"]) else 0
    greater_than = sum([allele_len >= cutoff for allele_len in allele_lens])
    
    return greater_than


def main(hits: pd.DataFrame, out_file: str, ensembl: str, constraint: str, omim: str, hpo: Optional[str] = None) -> None:
    # convert hits dataframe from long to wide format
    hits = pd.read_csv(hits)
    hits.rename({"case_trid": "trid"}, axis=1, inplace=True)
    hits_pivot = pivot_hits(hits)

    # make a column with maximum z score across samples
    z_score_cols = [col for col in hits_pivot.columns if "z_score" in col]
    for col in z_score_cols:
        hits_pivot[col] = [
            round(score, 3) if not pd.isnull(score) else None
            for score in hits_pivot[col]
        ]
    hits_pivot["max_z_score"] = hits_pivot[z_score_cols].max(axis=1)

    # filter out non-outliers
    print("Filter outliers")
    al_cols = [col for col in hits_pivot.columns if '_allele_len' in col]
    hits_pivot["outlier"] = hits_pivot.apply(lambda row: filter_outliers(row, al_cols), axis=1)
    hits_pivot = hits_pivot[hits_pivot["outlier"]]

    # make a column that sums the number of individuals carrying a particular repeat expansion
    hits_pivot['num_samples'] = hits_pivot.apply(lambda row: num_expanded(row, al_cols), axis=1)

    # convert hits to PyRanges object
    hits_pr = hits_to_pr(hits_pivot)

    # prep Ensembl gene GTF
    print("Prep Ensembl GTF")
    cols = ["gene_name", "gene_id", "gene_biotype", "Feature"]
    gene_gr = prepare_Ensembl_GTF(
        ensembl, cols=cols
    )

    # annotate hits with Ensembl genes
    print("Annotate against Ensembl genes")
    hits_gene = annotate_genes(hits_pr, pr.PyRanges(gene_gr))

    # annotate with gene constraint
    print("Add gnomAD gene constraint")
    constraint_cols = ["gene", "lof.oe_ci.upper", "lof.pLI"]
    constraint = pd.read_csv(constraint, sep="\t")[constraint_cols].dropna()
    hits_gene = add_constraint(constraint, hits_gene)

    # annotate with OMIM
    print("Add OMIM phenotype")
    omim = prepare_OMIM(f"{omim}/mim2gene.txt", f"{omim}/morbidmap.txt")
    hits_gene_omim = annotate_OMIM(hits_gene, omim)

    # annotate with HPO terms
    if hpo==None:
        print("No HPO terms supplied")
        hits_gene_omim["HPO"] = ""
    else:
        print("Add HPO terms")
        hpo = pd.read_csv(hpo, sep="\t")
        hits_gene_omim = add_hpo(hpo, hits_gene_omim)

    # group and aggregate gene columns
    hits_gene_omim = group_by_gene(hits_gene_omim)

    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        hits_gene_omim[col] = hits_gene_omim[col].apply(lambda genes: gene_set(genes))

    # shorten sample names (e.g. remove .m84090_240207_191948_s1.hifi_reads.bc2013.KL.GRCh38.aligned.haplotagged.trgt.sorted)
    if 'haplotagged' in al_cols[0]:
        new_al_cols = []
        for col in al_cols:
            new_col = col.split('.')[0]
            new_col = new_col + "_allele_len"
            new_al_cols.append(new_col)
            hits_gene_omim.rename({col: new_col}, inplace=True, axis=1)
        
        new_z_score_cols = []
        for col in z_score_cols:
            new_col = col.split('.')[0]
            new_col = new_col + "_z_score"
            new_z_score_cols.append(new_col)
            hits_gene_omim.rename({col: new_col}, inplace=True, axis=1)
    else:
        new_al_cols = al_cols
        new_z_score_cols = z_score_cols
    
    hits_gene_omim = hits_gene_omim[
        [
            "Chromosome",
            "Start",
            "End",
            "trid",
            "gene_name",
            "gene_id",
            "gene_biotype"]
            + constraint_cols
            + [
            "Feature",
            "control_range",
            "cutoff",
            "max_z_score",
            "num_samples"
        ]
        + new_al_cols
        + new_z_score_cols
        + ["OMIM_phenotype", "HPO"]
    ]
    hits_gene_omim = hits_gene_omim.rename(
        columns={"Feature": "feature", "Chromosome": "CHROM", "Start": "POS", "End": "END", "trid": "TRID", "gene": "gnomad_constraint_gene"}
    )

    # recode missing values
    hits_gene_omim["CHROM"] = hits_gene_omim["CHROM"].astype(str)
    hits_gene_omim.fillna(".", inplace=True)
    hits_gene_omim.replace({"-1": "."}, inplace=True)

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

    args = parser.parse_args()
    print("Annotating repeats...")
    main(args.repeats, args.output_file, args.ensembl_gtf, args.gnomad_constraint, args.OMIM_path, args.hpo)