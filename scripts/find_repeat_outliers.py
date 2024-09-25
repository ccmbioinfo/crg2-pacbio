import argparse
import pandas as pd
import numpy as np


def split_alleles(alleles: pd.Series) -> tuple:
    """
    Split alleles by comma and return lengths and sequences of shortest and longest alleles.
    Note that there may be more than two alleles called by TRGT at a locus in an individual (somaticism?)
    """
    alleles = alleles.split(",")
    # meth = meth.split(",")
    # mp = mp.split(",")
    allele_lens = [len(a) for a in alleles]
    short_allele_len, long_allele_len = min(allele_lens), max(allele_lens)
    short_index, long_index = allele_lens.index(min(allele_lens)), allele_lens.index(
        max(allele_lens)
    )
    short_allele, long_allele = alleles[short_index], alleles[long_index]
    # short_allele_meth, long_allele_meth = meth[short_index], meth[long_index]
    # short_allele_mp, long_allele_mp = mp[short_index], mp[long_index]

    # return short_allele, short_allele_len, short_allele_meth, short_allele_mp, long_allele, long_allele_len, long_allele_meth, long_allele_mp
    return (
        short_allele,
        short_allele_len,
        short_index,
        long_allele,
        long_allele_len,
        long_index,
    )


def split_meth_mp(
    meth_mp: pd.Series, short_index: pd.Series, long_index: pd.Series
) -> tuple:
    """
    Split methylation or motify purity and return value for shortest and longest allele
    """
    try:
        meth_mp_short = meth_mp.split(",")[short_index]
        meth_mp_long = meth_mp.split(",")[long_index]
    except AttributeError:
        meth_mp_short = np.nan
        meth_mp_long = np.nan
    return meth_mp_short, meth_mp_long


def calc_z_score(
    allele_type: str,
    allele_length: int,
    short_allele_len_mean: float,
    short_allele_len_std: float,
    long_allele_len_mean: float,
    long_allele_len_std: float,
    am: float,
    AM_mean: float,
    AM_std: float,
    mp: float,
    MP_mean: float,
    MP_std: float,
) -> tuple:
    """
    Calculate z-score for sample allele length, methylation, and motif purity compared to control distribution
    """

    means = {
        "long_allele_len_mean": long_allele_len_mean,
        "short_allele_len_mean": short_allele_len_mean,
        "AM_mean": AM_mean,
        "MP_mean": MP_mean,
    }
    stds = {
        "long_allele_len_std": long_allele_len_std,
        "short_allele_len_std": short_allele_len_std,
        "AM_std": AM_std,
        "MP_std": MP_std,
    }

    def calculate_z_score(value, mean, std):
        try:
            z_score = (value - mean) / std
        except ZeroDivisionError:
            z_score = (value - mean) / 0.0001  # handle division by zero
        return round(z_score, 3)

    z_score_len = calculate_z_score(
        allele_length, means[f"{allele_type}_len_mean"], stds[f"{allele_type}_len_std"]
    )
    try:
        z_score_am = calculate_z_score(
            float(am), means["AM_mean"] / 255, stds["AM_std"] / 255
        )  # CMH TRGT VCFs called with v0.4.0, AM a value between 0 and 255
    except:
        z_score_am = "."
    try:
        z_score_mp = calculate_z_score(float(mp), means["MP_mean"], stds["MP_std"])
    except:
        z_score_mp = "."

    return z_score_len, z_score_am, z_score_mp


def main(cases, dist, output_file):
    # get length of shortest and longest alleles
    # cases["short_allele"], cases["short_allele_len"], cases["short_allele_meth"], cases["short_allele_MP"], cases["long_allele"], cases["long_allele_len"], cases["long_allele_meth"], cases["long_allele_MP"], = zip(
    #     *cases.apply(lambda x: split_alleles(x.alleles, x.avg_methylation, x.motif_purity), axis=1)
    # )

    (
        cases["short_allele"],
        cases["short_allele_len"],
        cases["short_allele_index"],
        cases["long_allele"],
        cases["long_allele_len"],
        cases["long_allele_index"],
    ) = zip(*cases.apply(lambda x: split_alleles(x.alleles), axis=1))

    cases["short_allele_meth"], cases["long_allele_meth"] = zip(
        *cases.apply(
            lambda x: split_meth_mp(
                x.avg_methylation, x.short_allele_index, x.long_allele_index
            ),
            axis=1,
        )
    )

    cases["short_allele_MP"], cases["long_allele_MP"] = zip(
        *cases.apply(
            lambda x: split_meth_mp(
                x.motif_purity, x.short_allele_index, x.long_allele_index
            ),
            axis=1,
        )
    )

    # drop alleles, not needed for now?
    cases = cases.drop(columns=["alleles"]).copy()

    # convert dataframe from wide to long format
    # allele lengths
    cases_al = pd.melt(
        cases,
        id_vars=["trid", "sample"],
        value_vars=["short_allele_len", "long_allele_len"],
        var_name="allele_length_type",
        value_name="allele_length",
    )
    cases_al = cases_al.drop(columns=["trid", "sample"])

    # allele sequences
    cases_al_seq = pd.melt(
        cases,
        id_vars=["trid", "sample"],
        value_vars=["short_allele", "long_allele"],
        var_name="allele_type",
        value_name="allele_sequence",
    )
    cases_al_seq = cases_al_seq.drop(columns=["trid", "sample"])

    # methylation
    cases_meth = pd.melt(
        cases,
        id_vars=["trid", "sample"],
        value_vars=["short_allele_meth", "long_allele_meth"],
        var_name="allele_am_type",
        value_name="AM",
    )
    cases_meth = cases_meth.drop(columns=["trid", "sample"])

    # motif purity
    cases_mp = pd.melt(
        cases,
        id_vars=["trid", "sample"],
        value_vars=["short_allele_MP", "long_allele_MP"],
        var_name="allele_mp_type",
        value_name="MP",
    )

    cases_long = pd.concat([cases_al, cases_al_seq, cases_meth, cases_mp], axis=1)

    # merge case alleles to controls allele length distribution dataframe
    cases_long["case_trid"] = cases_long["trid"]
    cases_long["trid"] = cases_long["case_trid"].str.rsplit("_", n=1).str[0]
    dist["trid"] = dist["trid"].str.rsplit("_", n=1).str[0]
    merged = cases_long.merge(dist, on="trid", how="left")
    merged = merged.astype(
        {
            "short_allele_len_mean": float,
            "long_allele_len_mean": float,
            "short_allele_len_std": float,
            "long_allele_len_std": float,
            "cutoff_short": float,
            "cutoff_long": float,
            "MP": float,
            "AM_mean": float,
            "AM_std": float,
            "MP_mean": float,
            "MP_std": float,
        }
    )

    # calculate z_scores for each case allele compared to controls
    print("Calculating z scores")
    merged["z_score_len"], merged["z_score_AM"], merged["z_score_MP"] = zip(
        *merged.apply(
            lambda row: calc_z_score(
                row["allele_type"],
                row["allele_length"],
                row["short_allele_len_mean"],
                row["short_allele_len_std"],
                row["long_allele_len_mean"],
                row["long_allele_len_std"],
                row["AM"],
                row["AM_mean"],
                row["AM_std"],
                row["MP"],
                row["MP_mean"],
                row["MP_std"],
            ),
            axis=1,
        )
    )

    # split into  two dataframes: short alleles and long alleles
    merged_short = merged[merged["allele_type"] == "short_allele"][
        [
            "case_trid",
            "sample",
            "allele_type",
            "allele_length",
            "cutoff_short",
            "z_score_len",
            "range_short",
            "AM",
            "AM_mean",
            "AM_std",
            "z_score_AM",
            "MP",
            "MP_mean",
            "MP_std",
            "z_score_MP",
        ]
    ]
    merged_long = merged[merged["allele_type"] == "long_allele"][
        [
            "case_trid",
            "sample",
            "allele_type",
            "allele_length",
            "cutoff_long",
            "z_score_len",
            "range_long",
            "AM",
            "AM_mean",
            "AM_std",
            "z_score_AM",
            "MP",
            "MP_mean",
            "MP_std",
            "z_score_MP",
        ]
    ]

    # cat short and long alleles into one dataframe
    print("Concatentating short and long alleles")
    merged_short.rename(
        {"range_short": "control_range", "cutoff_short": "cutoff"}, inplace=True, axis=1
    )
    merged_long.rename(
        {"range_long": "control_range", "cutoff_long": "cutoff"}, inplace=True, axis=1
    )
    merged_cat = pd.concat([merged_short, merged_long], axis=0)
    
    # shorten sample names (e.g. remove .m84090_240207_191948_s1.hifi_reads.bc2013.KL.GRCh38.aligned.haplotagged.trgt.sorted)
    merged_cat["sample"] = merged_cat["sample"].apply(lambda x: x.split('.')[0])

    # remove unnecessary columns
    merged_cat = merged_cat[
        [
            "case_trid",
            "sample",
            "allele_type",
            "allele_length",
            "z_score_len",
            "control_range",
            "cutoff",
            "AM",
            "AM_mean",
            "AM_std",
            "z_score_AM",
            "MP",
            "MP_mean",
            "MP_std",
            "z_score_MP",
        ]
    ]
    merged_cat.rename({"allele_length": "allele_len"}, axis=1, inplace=True)

    # sort by z_score
    merged_cat = merged_cat.sort_values(by="z_score_len", ascending=False)

    # round
    merged_cat.round(2)

    # export
    merged_cat.to_csv(output_file, index=False)


if __name__ == "__main__":
    # if running from the command-line
    description = "Find outlier repeats in family vs controls"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--alleles_path",
        type=str,
        required=True,
        help="Repeat alleles for cases",
    )
    parser.add_argument(
        "--control_alleles",
        type=str,
        required=True,
        help="Control allele length distribution CSV",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )

    args = parser.parse_args()
    cases = pd.read_csv(
        args.alleles_path,
        sep=" ",
        compression="gzip",
        header=None,
        names=["trid", "sample", "motif_purity", "avg_methylation", "alleles"],
    )
    print("Loading control distributions")
    dist = pd.read_csv(args.control_alleles, sep="\t")
    print("Determining outlier repeats")
    main(cases, dist, args.output_file)
