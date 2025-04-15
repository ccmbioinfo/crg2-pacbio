import argparse
import pandas as pd
import pysam
import numpy as np

def resample_quantiles(counts, num_resamples=100):
    resamples = np.random.choice(counts, len(counts) * num_resamples).reshape(num_resamples, -1)
    return list(np.quantile(resamples, 0.95, axis=1))

def get_counts_with_finite_std(counts):
    return counts[:-1] + [counts[-1] + 0.1] if len(set(counts)) == 1 else counts

def get_cutoff(quantiles):
    return np.mean(quantiles) + max(1, np.std(quantiles))

def get_quantiles_cutoffs(allele_lens):
    quantiles = get_counts_with_finite_std(resample_quantiles(allele_lens))
    return get_cutoff(quantiles), [min(allele_lens), max(allele_lens)]

def process_variant(variant, samples):
    allele_dict = {
        "short_allele_len": [], "long_allele_len": [], "short_AM": [], "long_AM": [],
        "short_MP": [], "long_MP": []
    }
    for sample in samples:
        allele_lens = variant.samples[sample]['AL']
        AM = variant.samples[sample]['AM']
        try:
            MP = variant.samples[sample]["MP"]
        except KeyError:  # Fallback to AP for older TRGT version
            MP = variant.samples[sample]["AP"]
        if allele_lens[0] is not None:
            min_index, max_index = allele_lens.index(min(allele_lens)), allele_lens.index(max(allele_lens))
            allele_dict["short_allele_len"].append(min(allele_lens))
            allele_dict["long_allele_len"].append(max(allele_lens))
            allele_dict["short_AM"].append(AM[min_index])
            allele_dict["long_AM"].append(AM[max_index])
            allele_dict["short_MP"].append(MP[min_index])
            allele_dict["long_MP"].append(MP[max_index])
    return allele_dict

def calculate_stats(allele_dict):
    stats = {}
    for key in ["short_allele_len", "long_allele_len"]:
        stats[f"{key}_mean"] = np.mean(allele_dict[key])
        stats[f"{key}_std"] = np.std(allele_dict[key])
        stats[f"cutoff_{key.split('_')[0]}"], stats[f"range_{key.split('_')[0]}"] = get_quantiles_cutoffs(allele_dict[key])
        # get z-scores for all alleles
        if stats[f"{key}_std"] == 0:
            stats[f"{key}_std"] = 0.1
        stats[f"{key}_zscore"] = [round((x - stats[f"{key}_mean"]) / stats[f"{key}_std"], 2) for x in allele_dict[key]]

    # Calculate mean and std for average methylation and motif purity
    for metric in ["AM", "MP"]:
        combined = allele_dict[f"short_{metric}"] + allele_dict[f"long_{metric}"]
        combined = [x for x in combined if not pd.isnull(x)]
        stats[f"{metric}_mean"] = np.mean(combined) if combined else np.nan
        stats[f"{metric}_std"] = np.std(combined) if combined else np.nan
    
    return stats

def main(input_vcf, output_tsv):
    with pysam.VariantFile(input_vcf, 'r') as vcf_in:
        samples = list(vcf_in.header.samples)
        cutoff_dict = {
            "trid": [], "short_allele_len_mean": [], "long_allele_len_mean": [],
            "short_allele_len_std": [], "long_allele_len_std": [],
            "cutoff_short": [], "cutoff_long": [], "range_long": [], "range_short": [],
            "AM_mean": [], "AM_std": [], "MP_mean": [], "MP_std": [],
            "short_allele_len_zscore": [], "long_allele_len_zscore": [] 
        }

        for variant in vcf_in:
            trid = variant.info['TRID']
            print(f"Reading variant {trid}")
            allele_dict = process_variant(variant, samples)
            
            if not allele_dict["short_allele_len"]:  # Empty because genotype is missing for every sample
                stats = {key: np.nan for key in cutoff_dict.keys() if key != "trid"}
            else:
                stats = calculate_stats(allele_dict)
            
            cutoff_dict["trid"].append(trid)
            for key, value in stats.items():
                cutoff_dict[key].append(value)

        df = pd.DataFrame(cutoff_dict).round(2)
        df.to_csv(output_tsv, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TRGT multi-sample VCF and calculate allele statistics.")
    parser.add_argument("--input_vcf", help="Path to the input multi-sample TRGT VCF file")
    parser.add_argument("--output_tsv", help="Path to the output TSV file")
    args = parser.parse_args()

    main(args.input_vcf, args.output_tsv)
