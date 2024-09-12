# Calculating tandem repeat allele length distributions

The scripts in this folder generate per-chromosome TSV files that describe tandem repeat allele length, motif purity, and methylation ditributions.
The TSV will contain one tandem repeat per row. The columns are shown below. 

Note that I was using a conda environment to run this. You'll just need python3 with pandas and numpy installed.

Many of the python functions were originally written by Egor Dolzhenko from PacBio, which I have since adapted. See [this workflow](https://github.com/tandem-repeat-workflows/find-outlier-expansions/blob/main/find-outlier-expansions.ipynb) for original outlier detection code. 


| trid                       | short_allele_len_mean | short_allele_len_std | long_allele_len_mean | long_allele_len_std | MP_mean | MP_std | AM_mean | AM_std | cutoff_short | range_short | cutoff_long | range_long |
|----------------------------|-----------------------|----------------------|----------------------|---------------------|---------|--------|---------|--------|--------------|-------------|-------------|------------|
| chr22_10516072_10516291_TA | 230.28                | 8.03                 | 236.38               | 6.43                | 0.74    | 0.01   | 0.65    | 0.14   | 245.7        | [219, 247]  | 247.99      | [219, 249] |


1. generate_allele_db.sh: generates database of tandem repeat alleles, \<prefix\>.alleles.db.gz, from provided TRGT VCFs, and then splits the output by chromosome. The output TSVs contain one allele locus per sample per row. 
2. generate_allele_dist.sh: after running the script above, this script calculates  tandem repeat allele length, motif purity, and methylation ditributions per-chromosome for each locus and output a TSV file as shown above. The script calculates the allele length distributions for the short allele and long alleles separately, i.e. we aggregate the shortest allele for all samples and calculate sample statistics, and do the same for the longest allele for all samples. This may be useful in autosomal recessive TR expansion inheritance; a TR expansion may not be an outlier if the allele is present in the population, but if we consider the short and long alleles separately, the 'short' allele should be picked up as an outlier if there are few or no homozygous carriers in the population. The cutoffs (cutoff_short and cutoff_long) are calculated as the mean plus standard deviation of the 95th quantile for TR allele length (resampled 100 times). 