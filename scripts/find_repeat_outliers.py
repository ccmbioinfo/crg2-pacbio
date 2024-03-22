import argparse
import gzip
import itertools
import os
import pandas as pd
import numpy as np
from collections import namedtuple
from pathlib import Path

# based on Egor Dolzhenko's code in https://github.com/tandem-repeat-workflows/find-outlier-expansions/blob/main/find-outlier-expansions.ipynb

RepeatRec = namedtuple("RepeatRec", "sample short_allele long_allele")

def get_repeat_recs(path):
    def parse_alleles(group):
        alleles = list(line.decode("utf8").split() for line in group)
        alleles = [(rec[1], rec[2].split(",")) for rec in alleles]
        return alleles
    
    with gzip.open(path, "r") as file:
        for trid, group in itertools.groupby(file, key=lambda line: line.decode("utf8").split()[0]):
            alleles = parse_alleles(group)
            repeat_recs = [(s, [len(a) for a in als]) for s, als in alleles]
            repeat_recs = [RepeatRec(s, min(als), max(als)) for s, als in repeat_recs]
            
            yield trid, repeat_recs


def resample_quantiles(counts, num_resamples):
    """Based on https://github.com/Illumina/ExpansionHunterDenovo/blob/master/scripts/core/common.py"""
    resamples = np.random.choice(counts, len(counts) * num_resamples)
    resamples = np.split(resamples, num_resamples)

    resampled_quantiles = []
    for resample in resamples:
        quantile = np.quantile(resample, 0.95)
        resampled_quantiles.append(quantile)

    return resampled_quantiles


def get_counts_with_finite_std(counts):
    if len(set(counts)) == 1:
        return counts[:-1] + [counts[-1] + 0.1]
    return counts


def get_cutoff(quantiles):
    mean = np.mean(quantiles)
    std = max(1, np.std(quantiles))
    cutoff = mean + std
    return cutoff


def get_hits(allele_type, repeat_recs, case_ids):
    assert allele_type in ["long", "short"]
    allele_index = 2 if allele_type == "long" else 1
    cases, controls = {}, []
    for rec in repeat_recs:
        if rec.sample in case_ids:
            cases[rec.sample] = rec[allele_index]
        else:
            controls.append(rec[allele_index])
    if not len(controls) == 0:
        quantiles = resample_quantiles(controls, 100)
        quantiles = get_counts_with_finite_std(quantiles)
        cutoff = get_cutoff(quantiles)
    
        for case, allele_len in cases.items():
            mean_control = np.mean(controls)
            std_control = np.std(controls)
            z_score = (allele_len - mean_control) / std_control
            if allele_len > cutoff:
                yield case, allele_len, z_score, (min(controls), max(controls))
    else:
        # no controls carry this TRID
        for case, allele_len in cases.items():
            yield case, allele_len, np.nan, np.nan


def main(alleles_path, out_filepath, case_ids):    
    if not os.path.exists("repeat_outliers"):
        Path("repeat_outliers").mkdir(exist_ok=True) 
    alleles_path = Path(alleles_path).resolve(strict=True)
    HitRec = namedtuple("HitRec", "trid sample allele_type allele_len z_score control_range")
    hits = []
    for i, (trid, repeat_recs) in enumerate(get_repeat_recs(alleles_path)):    
        for hit_type in ["long", "short"]:
            for case, allele_len, z_score, control_range in get_hits(hit_type, repeat_recs, case_ids):
                hits.append(HitRec(trid, case, hit_type, allele_len, z_score, control_range))

    df = pd.DataFrame(hits)
    df.to_csv(out_filepath, index=False)


if __name__ == "__main__":
    # if running from the command-line
    description = "Find outlier repeats in family vs controls"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--alleles_path",
        type=str,
        required=True,
        help="Alleles from generate_allele_db",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )
    parser.add_argument(
        "--case_ids",
        type=str,
        required=True,
        help="Path to Ensembl gene GTF",
    )

    args = parser.parse_args()
    print("Determining outlier repeats")
    case_ids = args.case_ids
    case_ids = pd.read_table(case_ids)["sample"].to_list()
    main(args.alleles_path, args.output_file, case_ids)