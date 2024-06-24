import argparse
import os
import pandas as pd
import numpy as np


def resample_quantiles(counts: pd.Series, num_resamples: int) -> list:
    """Based on https://github.com/Illumina/ExpansionHunterDenovo/blob/master/scripts/core/common.py"""
    resamples = np.random.choice(counts, len(counts) * num_resamples)
    resamples = np.split(resamples, num_resamples)

    resampled_quantiles = []
    for resample in resamples:
        quantile = np.quantile(resample, 0.95)
        resampled_quantiles.append(quantile)

    return resampled_quantiles


def get_counts_with_finite_std(counts: pd.Series) -> pd.Series:
    if len(set(counts)) == 1:
        return counts[:-1] + [counts[-1] + 0.1]
    return counts


def get_cutoff(quantiles: list) -> str:
    mean = np.mean(quantiles)
    std = max(1, np.std(quantiles))
    cutoff = mean + std
    return cutoff


def split_alleles(alleles: pd.DataFrame) -> tuple:
    """
    Split alleles by comma and return lengths of shortest and longest alleles.
    Note that there may be more than two alleles called by TRGT at a locus in an individual (somaticism?)
    """
    allele_lens = [len(a) for a in alleles.split(",")]
    short_allele_len = min(allele_lens)
    long_allele_len = max(allele_lens)

    return short_allele_len, long_allele_len


def get_quantiles_cutoffs(trid: str, allele_type: str, alleles: pd.DataFrame) -> tuple:
    """
    Calculate the a length cutoff for a particular tandem repeat locus as defined by the mean plus std of 95% quantile after resampling
    """
    alleles_trid = alleles[alleles["trid"] == trid]
    min = alleles_trid[allele_type].min()
    max = alleles_trid[allele_type].max()
    quantiles = resample_quantiles(alleles_trid[allele_type], 100)
    quantiles = get_counts_with_finite_std(quantiles)
    cutoff_short = get_cutoff(quantiles)

    return cutoff_short, [min, max]

def main(alleles):
    print("Loading allele db")
    dirname = os.path.dirname(alleles)
    outfile = os.path.basename(alleles).replace(".db", ".distribution.tsv.gz")
    outpath = os.path.join(dirname, outfile)
    alleles = pd.read_csv(
        alleles,
        header=None,
        names=["trid", "sample", "motif_purity", "avg_methylation", "alleles"],
        sep=" ",
    )

    # get length of shortest and longest alleles
    alleles["short_allele_len"], alleles["long_allele_len"] = zip(
        *alleles["alleles"].map(lambda x: split_alleles(x))
    )

    try:
        # split motif purity column
        alleles[["MP1", "MP2"]] = alleles["motif_purity"].str.split(
            ",",
            expand=True,
        )
        # split average methylation column
        alleles[["AM1", "AM2"]] = alleles["avg_methylation"].str.split(
            ",",
            expand=True,
        )
    except:
        # on Y chromosome, MP is a single float value
        # for simplicity, just encode MP1/MP2 and AM1/AM2 as the same value
        alleles["MP1"] = alleles["motif_purity"]
        alleles["MP2"] = alleles["motif_purity"]
        alleles["AM1"] = alleles["avg_methylation"]
        alleles["AM2"] = alleles["avg_methylation"]

    alleles = alleles.replace(".", np.nan)
    alleles = alleles.astype({"MP1": float, "MP2": float})
    alleles = alleles.astype({"AM1": float, "AM2": float})

    # calculate mean motify purity and average methylation at locus
    # should we separate these out in the future (i.e. for short and long alleles)? do we need resampling?
    alleles["MP"] = alleles[["MP1", "MP2"]].mean(axis=1)
    alleles["AM"] = alleles[["AM1", "AM2"]].mean(axis=1)

    # group alleles by tandem repeat ID and extract statistics
    print("Grouping alleles")
    grouped_alleles = alleles.groupby("trid").agg(
        {
            "short_allele_len": ["mean", "std"],
            "long_allele_len": ["mean", "std"],
            "MP": ["mean", "std"],
            "AM": ["mean", "std"],
        }
    )
    grouped_alleles = grouped_alleles.reset_index()

    # calculate thresholds for downstream outlier processing
    print("Computing quantiles")
    grouped_alleles["cutoff_short"], grouped_alleles["range_short"] = zip(
        *grouped_alleles["trid"].map(
            lambda x: get_quantiles_cutoffs(x, "short_allele_len", alleles)
        )
    )
    grouped_alleles["cutoff_long"], grouped_alleles["range_long"] = zip(
        *grouped_alleles["trid"].map(
            lambda x: get_quantiles_cutoffs(x, "long_allele_len", alleles)
        )
    )

    grouped_alleles = grouped_alleles.round(2)

    columns = [
        "trid",
        "short_allele_len_mean",
        "short_allele_len_std",
        "long_allele_len_mean",
        "long_allele_len_std",
        "MP_mean",
        "MP_std",
        "AM_mean",
        "AM_std",
        "cutoff_short",
        "range_short",
        "cutoff_long",
        "range_long",
    ]

    grouped_alleles.to_csv(
        outpath, sep="\t", compression="gzip", header=columns, index=False
    )


if __name__ == "__main__":
    # if running from the command-line
    description = "Calculate tandem repeat allele length distributions"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--alleles_path",
        type=str,
        required=True,
        help="Alleles from generate_allele_db.py",
    )

    args = parser.parse_args()
    main(args.alleles_path)
