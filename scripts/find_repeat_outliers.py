import argparse
import pandas as pd
import numpy as np
import pysam


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

def get_rank(zscores, zscore):
    # get the rank of a zscore in a list
    if zscore == np.nan:
        return np.nan
    else:
        sorted_list = sorted(zscores, reverse=True)
        # Find the rightmost/largest index of zscore in the sorted list
        for i in range(len(sorted_list)-1, -1, -1):
            if sorted_list[i] == zscore:
                return i + 1
        return None

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

def sample_vcf_to_dict(case_vcf, LPS_dict):
    # first, make a dictionary of the case(s) alleles (sample_dict[sample][trid] = [short_allele_len, long_allele_len, short_allele_lps, long_allele_lps, AM_short_allele, AM_long_allele, MP_short_allele, MP_long_allele])
    sample_dict = {}
    with pysam.VariantFile(case_vcf, 'r') as vcf_in:
        samples = vcf_in.header.samples
        for variant in  vcf_in:
            trid = variant.info['TRID']
            if len(trid.split("_")) == 4: # older versions of TRGT include motif in TRID
                trid = variant.info['TRID'].rsplit("_", 1)[0]
            for sample in samples: # may be multiple cases/family members
                sample_prefix = sample.split(".")[0]
                # extract allele lengths
                allele_lens = variant.samples[sample]['AL']
                min_index, max_index = allele_lens.index(min(allele_lens)), allele_lens.index(max(allele_lens))
                # extract allele LPS
                short_allele_len, long_allele_len = allele_lens[min_index], allele_lens[max_index]
                lps = LPS_dict[sample_prefix].get(trid, None)
                if lps: 
                    try:
                        short_allele_lps = LPS_dict[sample_prefix][trid][min_index + 1]
                        long_allele_lps = LPS_dict[sample_prefix][trid][max_index + 1]
                    except KeyError:
                        try:
                            short_allele_lps = LPS_dict[sample_prefix][trid][max_index + 1]
                            long_allele_lps = LPS_dict[sample_prefix][trid][max_index + 1]
                        except KeyError: # for some loci, two different alleles of same length, and LPS output only lists allele 2 
                            try:
                                short_allele_lps = LPS_dict[sample_prefix][trid][2]
                                long_allele_lps = LPS_dict[sample_prefix][trid][2]
                            except KeyError:
                                try:
                                    short_allele_lps, long_allele_lps = LPS_dict[sample_prefix][trid][1], LPS_dict[sample_prefix][trid][1]
                                except KeyError:
                                    short_allele_lps, long_allele_lps = np.nan, np.nan
                else: 
                    short_allele_lps, long_allele_lps = np.nan, np.nan
                # extract allele methylation
                AM = variant.samples[sample]['AM']
                try:
                    AM = [round(am, 2) for am in AM]
                    AM_short_allele, AM_long_allele = AM[min_index], AM[max_index]
                except TypeError: # methylation missing
                    AM_short_allele, AM_long_allele = np.nan, np.nan
                # extract motif purity
                try:
                    MP = variant.samples[sample]["MP"]
                except KeyError:  # fallback to AP for older TRGT version
                    MP = variant.samples[sample]["AP"]
                
                try:
                    MP = [round(mp, 2) for mp in MP]
                    MP_short_allele, MP_long_allele = MP[min_index], MP[max_index]
                except TypeError: # methylation missing
                    MP_short_allele, MP_long_allele = np.nan, np.nan


                if sample_prefix not in sample_dict:
                    sample_dict[sample_prefix] = {}

                sample_dict[sample_prefix][trid] = [short_allele_len, long_allele_len, short_allele_lps, long_allele_lps, AM_short_allele, AM_long_allele, MP_short_allele, MP_long_allele]
    
    return sample_dict



# iterate through background variants and compute allele stats
# for each TRID, calculate stats
# add TRID stats to sample dictionary (to output to df)
def main(case_vcf, case_lps, control_vcf, output_file):
    LPS_dict = {}
    case_lps = pd.read_csv(case_lps, sep="\t")
    for _, row in case_lps.iterrows(): 
        # create dictionary in the format {sample: {TRID: {allele1: lps_len, allele2: lps_len}}} 
        sample = row["sample"]
        trid = row.trid
        if len(trid.split("_")) == 4: 
            trid = trid.rsplit("_", 1)[0]
        allele = row.allele
        lps = row.lps_len
        if sample not in LPS_dict: 
            LPS_dict[sample] = {}
        if trid not in LPS_dict[sample]: 
            LPS_dict[sample][trid] = {allele: lps}
        else:
            LPS_dict[sample][trid][allele] = lps 
    sample_dict = sample_vcf_to_dict(case_vcf, LPS_dict)
    with pysam.VariantFile(control_vcf, 'r') as vcf_in:
        samples = list(vcf_in.header.samples)
        sample_stat_dict = {
            "trid": [],
            "sample": [],
            "allele_type": [],
            "allele_len": [], 
            "z_score_len": [],
            "z_score_len_rank": [], 
            "allele_len_mean": [],
            "allele_len_std": [], 
            "cutoff": [], 
            "range": [], 
            "distance_to_cutoff": [],
            "AM": [],
            "AM_mean": [],
            "AM_std": [],
            "z_score_AM": [],
            "MP": [],
            "MP_mean": [],
            "MP_std": [],
            "z_score_MP": [],
            "LPS": [],
        }

        for variant in  vcf_in:
            trid = variant.info['TRID']
            if len(trid.split("_")) == 4: # older versions of TRGT include motif in TRID
                trid = variant.info['TRID'].rsplit("_", 1)[0]
            #print(f"Reading variant {trid}")
            allele_dict = process_variant(variant, samples) # extract allele lengths, methylation, and motif purity for all samples at this locus
            
            try:
                stats = calculate_stats(allele_dict) # calculate allele length, methylation, and motif purity distributions
            except IndexError:   # genotype is missing for every sample
                stats_keys = 'short_allele_len_mean', 'short_allele_len_std', 'cutoff_short', 'range_short', 'short_allele_len_zscore', 'long_allele_len_mean', 'long_allele_len_std', 'cutoff_long', 'range_long', 'long_allele_len_zscore', 'AM_mean', 'AM_std', 'MP_mean', 'MP_std'
                stats = {key: np.nan for key in stats_keys if key not in ["trid"]}
                
            for case in sample_dict: # pull repeat stats from controls and calculate sample-specific allele length Z scores
                for allele in ["short", "long"]:
                    if trid not in sample_dict[case]:
                        continue
                    sample_stat_dict["trid"].append(trid)
                    sample_stat_dict["sample"].append(case)
                    sample_stat_dict["allele_type"].append(allele)
                    allele_len_mean = stats[f"{allele}_allele_len_mean"]
                    sample_stat_dict["allele_len_mean"].append(allele_len_mean)
                    allele_len_std = stats[f"{allele}_allele_len_std"]
                    sample_stat_dict["allele_len_std"].append(allele_len_std)
                    cutoff = stats[f"cutoff_{allele}"]
                    sample_stat_dict["cutoff"].append(cutoff)
                    sample_stat_dict["range"].append(stats[f"range_{allele}"])
                    
                    if allele == "short":
                        allele_len = sample_dict[case][trid][0]
                        LPS = sample_dict[case][trid][2]
                        AM =  sample_dict[case][trid][4]
                        MP = sample_dict[case][trid][6]
                    else:
                        allele_len = sample_dict[case][trid][1]
                        LPS = sample_dict[case][trid][3]
                        AM =  sample_dict[case][trid][5]
                        MP = sample_dict[case][trid][7]

                    sample_stat_dict["allele_len"].append(allele_len)
                    sample_stat_dict["LPS"].append(LPS)
                    try:
                        zscore_len = (allele_len - allele_len_mean) / allele_len_std
                        dist_to_cutoff = allele_len - cutoff
                    except TypeError: # no genotype information available in this sample for this allele
                        zscore_len = np.nan
                        dist_to_cutoff = np.nan
                    except ZeroDivisionError: # allele length std is 0
                        zscore_len = (allele_len - allele_len_mean) / 0.1
                        dist_to_cutoff = allele_len - cutoff
                        
                    sample_stat_dict["z_score_len"].append(zscore_len)
                    control_z_scores = stats[f"{allele}_allele_len_zscore"]
                    try:
                        control_z_scores.append(zscore_len)
                        z_score_len_rank = get_rank(control_z_scores, zscore_len)
                    except AttributeError: # no control z-scores available
                        z_score_len_rank = np.nan
                    sample_stat_dict["z_score_len_rank"].append(z_score_len_rank)  
                    sample_stat_dict["distance_to_cutoff"].append(dist_to_cutoff)  
                    # add methylation and motif purity information 
                    sample_stat_dict["AM"].append(AM)
                    sample_stat_dict["AM_mean"].append(stats["AM_mean"])
                    sample_stat_dict["AM_std"].append(stats["AM_std"])
                    try:
                        zscore_AM = (AM - stats["AM_mean"]) / stats["AM_std"]
                    except: 
                        zscore_AM =  (AM - stats["AM_mean"]) / 0.1
                    sample_stat_dict["z_score_AM"].append(zscore_AM)
                    sample_stat_dict["MP"].append(MP)
                    sample_stat_dict["MP_mean"].append(stats["MP_mean"])
                    sample_stat_dict["MP_std"].append(stats["MP_std"])
                    try:
                        zscore_MP = (MP - stats["MP_mean"]) / stats["MP_std"]
                    except: 
                        zscore_MP =  (MP - stats["MP_mean"]) / 0.1
                    sample_stat_dict["z_score_MP"].append(zscore_MP)
        sample_stat_df = pd.DataFrame.from_dict(sample_stat_dict).round(2)
        sample_stat_df.to_csv(output_file, sep="\t", index=False)
        return sample_stat_df





if __name__ == "__main__":
    # if running from the command-line
    description = "Find outlier repeats in family vs controls"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--case_vcf",
        type=str,
        required=True,
        help="Case multi-sample TRGT VCF",
    )
    parser.add_argument(
        "--case_lps",
        type=str,
        required=True,
        help="LPS values calculated by TRGT-LPS",
    )
    parser.add_argument(
        "--control_vcf",
        type=str,
        required=True,
        help="Control multi-sample TRGT VCF",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )

    args = parser.parse_args()
    print("Determining outlier repeats")
    main(args.case_vcf, args.case_lps,  args.control_vcf, args.output_file)