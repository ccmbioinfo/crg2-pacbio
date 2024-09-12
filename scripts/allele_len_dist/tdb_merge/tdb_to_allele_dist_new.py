import numpy as np
import pandas as pd
import tdb

def get_observed_allele_lengths(data, samples):
    """
    Given a tdb and a list of sample names
    return a DataFrame of LocusID, allele_number, allele_length, sample_name
    """
    observed_allele_lengths = []
    for name in samples:
        print(name)
        sample_table = data["sample"][name]
        # Create an index of the alleles present in the sample
        idx_table = pd.MultiIndex.from_frame(sample_table[["LocusID", "allele_number"]])
        # Subset the alleles and pull their lengths
        result = alleles.loc[idx_table, ["allele_length"]].drop_duplicates()
        print("Index and allele lens", len(idx_table), len(result))
        result['sample_name'] = name
        avg_meth = [am for am in sample_table["average_methylation"].values]
        print(len(avg_meth))
        result['average_methylation'] = avg_meth
        observed_allele_lengths.append(result)

    # Sorting and indexing here speeds up calculating the distribution per-locus later
    return (pd.concat(observed_allele_lengths)
                .reset_index()
                .sort_values("LocusID")
                .set_index("LocusID"))

def resample_quantiles(counts: pd.Series, num_resamples: int) -> list:
    """Based on https://github.com/Illumina/ExpansionHunterDenovo/blob/master/scripts/core/common.py"""
    resamples = np.random.choice(counts, len(counts) * num_resamples)
    resamples = np.split(resamples, num_resamples)
    resampled_quantiles = np.quantile(resamples, 0.95, axis=1)
    mean_quantiles = np.mean(resampled_quantiles)
    std_quantiles = np.std(resampled_quantiles)

    return mean_quantiles, std_quantiles

def round_to_sig_figs(x, sig_figs=3):
    """Round float to specified sig figs, default is three"""
    if x == 0:
        return 0
    elif pd.isna(x):
        return np.nan
    return round(x, sig_figs - int(np.floor(np.log10(abs(x)))) - 1)

#data = tdb.load_tdb("/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/CMH_TRGT/tdbs/final_tdb/final.tdb")
data = tdb.load_tdb("/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/data/TCAG_TRGT/tdbs/level2_tdbs/final_merge.tdb")
samples = [_ for _ in data["sample"].keys()]

# Only need to do this once so we can quickly grab the subset of alleles for each sample table
alleles = data["allele"].set_index(["LocusID", "allele_number"])
allele_lengths = get_observed_allele_lengths(data, samples)

# Some loci only have a single observed allele (homozygous sites)
# We can `drop_duplicates` to remove the redundant information
allele_lengths = (allele_lengths
                           .reset_index()
                           .drop_duplicates()
                           .set_index("LocusID"))

# Quick little summary stats
print(f"Observed {len(allele_lengths):,} alleles in {len(samples):,} controls over {allele_lengths.index.nunique():,} loci")

# Calculate statistics for allele lengths and methylation
grouped = allele_lengths.groupby("LocusID")
al_stats = grouped["allele_length"].agg(
    AL_95=lambda x: resample_quantiles(x, 10),
    AL_range=lambda x: [x.min(), x.max()]
)
al_stats = al_stats.join(pd.DataFrame({
    'AL_95_mean': al_stats['AL_95'].apply(lambda x: x[0]),
    'AL_95_std': al_stats['AL_95'].apply(lambda x: x[1])
})).drop('AL_95', axis=1)
am_stats = grouped["average_methylation"].agg(['mean', 'std']).rename(columns={'mean': "AM_mean", 'std': "AM_std"})

# Combine statistics and merge with locus data
allele_lengths_loci = (
    data["locus"]
    .merge(pd.concat([al_stats, am_stats], axis=1), on="LocusID")
    .assign(
        AL_cutoff=lambda df: df["AL_95_mean"] + df["AL_95_std"],
        trid=lambda df: df.apply(lambda x: f"{x['chrom']}_{x['start']}_{x['end']}", axis=1)
    )
    [["LocusID", "trid", "AL_range", "AL_95_mean", "AL_95_std", "AM_mean", "AM_std"]]
)

# Round float columns to three significant figures
float_columns = ["AL_95_mean", "AL_95_std", "AM_mean", "AM_std"]
allele_lengths_loci[float_columns] = allele_lengths_loci[float_columns].applymap(lambda x: round_to_sig_figs(x, 3))

# Export results
allele_lengths_loci.to_csv("test.csv", index=False)