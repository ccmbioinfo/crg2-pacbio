import argparse
import pandas as pd
import numpy as np


def split_alleles(alleles: pd.DataFrame) -> tuple:
    """
    Split alleles by comma and return lengths and sequences of shortest and longest alleles.
    Note that there may be more than two alleles called by TRGT at a locus in an individual (somaticism?)
    """
    alleles = alleles.split(",")
    allele_lens = [len(a) for a in alleles]
    short_allele_len = min(allele_lens)
    long_allele_len = max(allele_lens)
    short_index = alleles.index(min(alleles, key=len))
    long_index = alleles.index(max(alleles, key=len))
    short_allele = alleles[short_index]
    long_allele = alleles[long_index]

    return short_allele, short_allele_len, long_allele, long_allele_len


def calc_z_score(allele_type: str, allele_length: int, short_allele_len_mean: float, short_allele_len_std: float, long_allele_len_mean: float, long_allele_len_std: float) -> float:
    """
    Calculate z-score for sample allele length compared to control distribution
    """
    if allele_type == "short_allele":
        if short_allele_len_std == 0:
            std = 0.0001
        else:
            std = short_allele_len_std
        allele_len_mean = short_allele_len_mean
    else:
        if long_allele_len_std == 0:
            std = 0.0001
        else:
            std = long_allele_len_std
        allele_len_mean = long_allele_len_mean
    z_score = (allele_length - allele_len_mean)/std

    return z_score

def main(cases, dist, output_file):
    # get length of shortest and longest alleles
    cases["short_allele"], cases["short_allele_len"], cases["long_allele"], cases["long_allele_len"] = zip(
        *cases["alleles"].map(lambda x: split_alleles(x))
    )

    # drop motif purity and methylation, not needed for now?
    cases = cases.drop(columns=["motif_purity", "avg_methylation", "alleles"]).copy()

    # convert dataframe from wide to long format 
    # allele lengths
    cases_al = pd.melt(cases, id_vars=['trid', 'sample'], value_vars=['short_allele_len', 'long_allele_len'], var_name="allele_length_type", value_name="allele_length")
    cases_al = cases_al.drop(columns=["trid", "sample"])
    # allele sequences
    cases_al_seq =  pd.melt(cases, id_vars=['trid', 'sample'], value_vars=['short_allele', 'long_allele'], var_name="allele_type", value_name="allele_sequence")
    cases_long = pd.concat([cases_al, cases_al_seq], axis=1)

    # merge case alleles to controls allele length distribution dataframe
    cases_long["case_trid"] = cases_long["trid"]
    cases_long["trid"] = cases_long["case_trid"].apply(lambda x: x.rsplit('_', 1)[0])
    dist["trid"] = dist["trid"].apply(lambda x: x.rsplit('_', 1)[0])
    merged = cases_long.merge(dist, on="trid", how="left")
    merged = merged.astype({"short_allele_len_mean": float,
                                            "long_allele_len_mean": float,
                                            "short_allele_len_std": float,
                                            "long_allele_len_std": float,
                                            "cutoff_short": float,
                                            "cutoff_long": float,
                                            })
    print(len(merged))
    # calculate z_scores for each case allele compared to controls
    print("Calculating z scores")
    merged["z_score"] = merged.apply(lambda row: calc_z_score(row["allele_type"],
                                                            row["allele_length"],
                                                            row["short_allele_len_mean"],
                                                            row["short_allele_len_std"],
                                                            row["long_allele_len_mean"],
                                                            row["long_allele_len_std"]), axis=1)
    
    # split into  two dataframes: short alleles and long alleles
    merged_short = merged[merged["allele_type"] == "short_allele"][["case_trid", "sample", "allele_type", "allele_length", "cutoff_short", "z_score", "range_short"]]
    merged_long = merged[merged["allele_type"] == "long_allele"][["case_trid", "sample", "allele_type", "allele_length", "cutoff_long", "z_score", "range_long"]]

    # cat short and long alleles into one dataframe
    merged_short.rename({'range_short': "control_range", "cutoff_short": "cutoff"}, inplace=True, axis=1)
    merged_long.rename({'range_long': "control_range", "cutoff_long": "cutoff"}, inplace=True, axis=1)
    merged_cat = pd.concat([merged_short, merged_long], axis=0)

    # remove unnecessary columns
    merged_cat = merged_cat[["case_trid", "sample", "allele_type", "allele_length", "z_score", "control_range", "cutoff"]]
    merged_cat.rename({"allele_length": "allele_len"}, axis=1, inplace=True)
    
    # sort by z_score
    merged_cat = merged_cat.sort_values(by="z_score", ascending=False)

    # replace infinite z scores (std of zero in controls) with fixed value to facilitate filtering
    merged_cat.replace({np.inf: 10}, inplace=True)

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
    cases = pd.read_csv(args.alleles_path, sep=" ", compression="gzip",
                    header=None, names=["trid", "sample", "motif_purity", "avg_methylation", "alleles"])
    dist = pd.read_csv(args.control_alleles, sep="\t")
    print("Determining outlier repeats")
    main(cases, dist, args.output_file)
