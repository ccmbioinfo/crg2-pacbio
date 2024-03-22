import argparse
from datetime import date
import numpy as np
import pandas as pd
import pyranges as pr

def pivot_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert hits dataframe from long to wide format so each sample is associated with two columns:
    [sample]_allele_len, and [sample]_z_score
    """
    hits = df.copy()
    hits["trid"] = hits["trid"] + "_" + hits["allele_type"]
    hits = hits.drop("allele_type", axis=1)
    hit_pivot = hits.pivot(
        index=["trid", "control_range"],
        columns="sample",
        values=["allele_len", "z_score"],
    ).reset_index()

    hit_pivot.columns = (
        hit_pivot.columns.get_level_values(1)
        + "_"
        + hit_pivot.columns.get_level_values(0)
    )
    hit_pivot = hit_pivot.rename(
        columns={"_trid": "trid", "_control_range": "control_range"}
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

def add_hpo(hpo: pd.DataFrame, gene: str) -> list:
    """
    Add gene-based HPO terms
    """
    try:
        genes = gene.split(";")
    except AttributeError:
        return np.nan
    terms = []
    for gene in genes:
        try:
            term = str(hpo[hpo["ensembl_gene_id"] == gene]["Features"].values[0])
            term = term.replace("; ", ";").split(";")
            term = list(set(term))
            for t in term:
                terms.append(t)
        except IndexError:
            pass
    if len(terms) == 0:
        return np.nan
    else:
        terms = ",".join(terms)
        return terms


def main(hits: pd.DataFrame, out_file: str, ensembl: str, constraint: str, omim: str, hpo: str) -> None:
    # convert hits dataframe from long to wide format
    hits = pd.read_csv(hits)
    hits_pivot = pivot_hits(hits)
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

    # group and aggregate gene columns
    hits_gene = group_by_gene(hits_gene)

    # annotate with OMIM
    print("Add OMIM phenotype")
    omim = prepare_OMIM(f"{omim}/mim2gene.txt", f"{omim}/morbidmap.txt")
    hits_gene_omim = annotate_OMIM(hits_gene, omim)

    # annotate with HPO terms
    print("Add HPO terms")
    hpo = pd.read_csv(hpo, sep="\t")
    hits_gene_omim["HPO"] = [
            add_hpo(hpo, gene) for gene in hits_gene_omim["gene_id"].values
        ]

    # make a column with maximum z score across samples
    print("Format dataframe")
    z_score_cols = [col for col in hits_gene_omim.columns if "z_score" in col]
    for col in z_score_cols:
        hits_gene_omim[col] = [
            round(score, 3) if not pd.isnull(score) else None
            for score in hits_gene_omim[col]
        ]
    hits_gene_omim["max_z_score"] = hits_gene_omim[z_score_cols].max(axis=1)

    # make a column that sums the number of individuals carrying a particular repeat expansion
    al_cols = [col for col in hits_gene_omim.columns if '_allele_len' in col]
    hits_gene_omim['num_samples'] = hits_gene_omim[al_cols].notna().sum(axis=1)

    # column cleanup
    for col in ["gene_name", "gene_id", "gene_biotype", "Feature"]:
        hits_gene_omim[col] = hits_gene_omim[col].apply(lambda genes: gene_set(genes))

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
            "max_z_score",
            "num_samples"
        ]
        + al_cols
        + z_score_cols
        + ["OMIM_phenotype", "HPO"]
    ]
    hits_gene_omim = hits_gene_omim.rename(
        columns={"Feature": "feature", "Chromosome": "chr", "Start": "start", "End": "end", "gene": "gnomad_constraint_gene"}
    )

    # recode missing values
    hits_gene_omim["chr"] = hits_gene_omim["chr"].astype(str)
    hits_gene_omim.fillna(".", inplace=True)
    hits_gene_omim.replace({"-1": "."}, inplace=True)

    # write to file 
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    out_file = out_file.replace(".csv", "")
    hits_gene_omim.to_csv(f"{out_file}.csv", index=False)
    hits_gene_omim.to_csv(f"{out_file}_{today}.csv", index=False)



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
        required=True,
        help="Path to HPO terms file",
    )

    args = parser.parse_args()
    print("Annotating repeats...")
    main(args.repeats, args.output_file, args.ensembl_gtf, args.gnomad_constraint, args.OMIM_path, args.hpo)