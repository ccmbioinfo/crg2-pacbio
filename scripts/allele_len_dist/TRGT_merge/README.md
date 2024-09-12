# TRGT merge and statistics scripts

This folder contains scripts for merging and analyzing TRGT (Tandem Repeat Genotyping Tool) output files. Below is a description of each script, its purpose, and its arguments.

## 1. TRGT_merge_VCFs.sh

This Bash script merges multiple TRGT VCF files into a single multi-sample VCF file using `trgt merge`. Note that the TRGT version must be 1.1.0 or greater (prior versions do not include this functionality).

### Usage:
`sbatch TRGT_merge_VCFs.sh /path/to/TRGTv1.1.0/trgt  /path/to/VCFs/ /path/to/reference/human_GRCh38_no_alt_analysis_set.fasta multi_sample_TRGT`

## 2. parse_multi-sample_TRGT_VCF.sh
This python script takes a multi-sample TRGT VCF and calculates allele length, methylation, and motif purity statistics at each tandem repeat locus and outputs a TSV file as shown above. The script calculates the allele length distributions for the short allele and long alleles separately, i.e. we aggregate the shortest allele for all samples and calculate sample statistics, and do the same for the longest allele for all samples. This may be useful in autosomal recessive TR expansion inheritance; a TR expansion may not be an outlier if the allele is present in the population, but if we consider the short and long alleles separately, the 'short' allele should be picked up as an outlier if there are few or no homozygous carriers in the population. The cutoffs (cutoff_short and cutoff_long) are calculated as the mean plus standard deviation of the 95th quantile for TR allele length (resampled 100 times).

| trid                       | short_allele_len_mean | long_allele_len_mean | short_allele_len_std | long_allele_len_std | cutoff_short | cutoff_long | range_long | range_short | AM_mean | AM_std | MP_mean | MP_std |
|----------------------------|-----------------------|----------------------|----------------------|---------------------|---------|--------|---------|--------|--------------|-------------|-------------|------------|
| chr22_10516072_10516291| 230.23                 | 234.83                 | 8.69               | 8.76            | 244.71    | 249.13   | [217, 258]   | [216, 258]   | 0.51       | 0.27   |  0.74      | 0.01|

### Usage: 
`sbatch parse_multi-sample_TRGT_VCF.sh --input_vcf TRGT_multi-sample.vcf.gz --output_tsv TRGT_allele_stats.tsv`
