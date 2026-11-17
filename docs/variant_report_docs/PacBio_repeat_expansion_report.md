# Tandem repeat expansion report

Madeline Couse, Giulia Del Gobbo, Christian Marshall, Kristin Kernohan

**Version 2025-10**

The tandem repeat (TR) expansion workflow (original version) was designed by Egor Dolzhenko and Tom Mokveld and adapted by Madeline Couse and Giulia Del Gobbo. Briefly, TRGT (v2.1.0) is used to genotype TRs from the [adotto repeat catalog](https://zenodo.org/record/7987365#.ZHY9TOzMJAc) (937,122 repeat loci) for each sample, and TR lengths are compared to control TR allele length distributions to identify outlying expansions.

**Control tandem repeat allele length distributions**

We calculated the TR allele length, average methylation, and motif purity distributions for the adotto repeat catalog from 1,436 individuals from the Children’s Mercy GA4K cohort (note – this cohort does have affected individuals). We calculated the distributions for the short allele and long alleles separately, i.e. we aggregated the shortest allele for all samples and calculated sample statistics and did the same for the longest allele for all samples. This may be useful in autosomal recessive TR expansion inheritance; a TR expansion may not be an outlier if the allele is present in the population, but if we consider the short and long alleles separately, the 'short' allele should be picked up as an outlier if there are few or no homozygous carriers in the population. The allele length cutoffs for the short and long alleles, used to help determine repeat alleles expanded relative to the control population, are calculated as the mean plus standard deviation of the 95th quantile for TR allele length (resampled 100 times).

**Finding TR expansions**

We define tandem repeat expansions as tandem repeat alleles whose length is greater than the cutoffs defined above and has a Z score of greater than or equal to three. The report only contains tandem repeat loci where at least one family member possesses an outlier expansion. Note that there may be a small number of loci that were missing from the control cohort, where the cutoff is annotated as ‘.’.

**Longest pure segment (LPS) length**

LPS length is the longest stretch of uninterrupted repeat motifs at a TR locus, for example, the longest continuous sequence of GAA repeats. Pathogenic repeats exhibit higher LPS length variance compared to the rest of repeat catalog. Interruptions in repeats, which reduce the length of the LPS, tend to stabilize alleles.

**Validation**

We validated this pipeline on expanded alleles at 11 disease-associated loci in Coriell samples. We found that sorting the outlier report by LPS length ranked 6 of 11 pathogenic alleles (*FXN, DMPK, FMR1, HTT, ATN1,* and *ATXN3*) in the top 10 repeat expansions (Table 1). *AR* and *ATXN1* alleles were ranked in the top 20 after several filters were applied. The *PABPN1* repeat expansion was challenging, ranking 1,816: \>=12 GCN repeats are pathogenic in a dominant manner, while 11 repeats are pathogenic in a recessive manner. The pathogenic allele LPS was 11 for the long allele, or 11 repeat copies. The relatively short length of pathogenic alleles at this locus means that prioritization by LPS will not be effective. Two alleles were not called outliers and therefore did not make it into the report: one was a relatively small expansion in *C9ORF72* which may not be pathogenic based on the literature, and the other was an expansion of a non-reference motif in *RFC1*. Pathogenic expansion of non-reference motifs at loci where expanded reference motifs are benign will be challenging to detect with this pipeline since the adotto catalog is reference-based.

| **Gene** | **Repeat outlier report rank when sorted by “max_LPS” col, unfiltered** | **Repeat outlier report rank when sorted by “max_LPS”, filtered** |
| --- | --- | --- |
| FXN | 2 | 1 |
| DMPK | 3 | 2 |
| FMR1 | 4 | 3 |
| HTT | 6 | 5 |
| ATN1 | 7 | 7 |
| ATXN3 | 12 | 9 |
| AR | 36 | 14 |
| ATXN1 | 77 | 24 |
| PABPN1 | 1816 | 151 |

**Table 1.** Repeat outlier report rank when outliers are sorted by “max_LPS” column, with and without filtering. Total number of outliers was 3,670. Filters applied: gene_id != “.”, C4R_outlier_count \<= 5, cmh003409_01_LPS_rank \<= 5. In a family-based report, we can filter down even further based on the number of samples carrying a repeat expansion. **Note that all pathogenic alleles here were spiked into one random sample, and the report rank reflects this – for example, if the *ATXN3* allele was the only spiked-in allele, it would have ranked 7<sup>th</sup> in the unfiltered report, and 4<sup>th</sup> in the filtered report**.

**TR constraint**

TR constraint measures the observed-to-expected ratio for length of the longest pure segment at a TR locus. Pathogenic repeats exhibit higher LPS length variance than expected compared to rest of repeat catalog (TR_constraint \> 1, though protein-coding loci are under some constraint, TR_constraint \< 1), therefore novel disease-causing loci are expected be highly polymorphic in general population (Danzi et al 2025).

**Analysis strategy**

The report is sorted by the ‘max_LPS’ column, which corresponds to the maximum LPS length across all family members for a given allele. If you have multiple affected family members, consider setting **num_outlier_samples** equal to the number of affected family members. Note that allelic dropout may occur for expanded alleles, or genotypes may be missing due to quality issues at a locus.

Recommended filters:

- C4R_outlier_count ≤ 5
- gene_id != “.”
- \<proband\>\_LPS_rank ≤ 10 (could reduce down to 5)
- gene_biotype includes “protein_coding” (NB: will not capture lncRNA associated with repeat expansions, such as ATXN8OS – perhaps best to start without this filter)
- For known disease genes: HPO not equal to “.”

Outliers are already sorted by max_LPS, so pay particular attention to the top outliers after filtering. Things to consider when interpreting outliers to identify candidates:

- LPS rank and z-score of LPS
- Motif: homopolymer expansions aren’t (yet) known to cause STR diseases, expansions of certain motifs like intronic GAAn, coding CAGn or GCNn contribute to numerous STR diseases
- Gene context
- If present in C4R, check the lengths compared to your test sample and whether other C4R are/aren’t affected with similar presentation (though note that this length outlier method doesn’t consider if there are different motifs in the population)
- Gene (especially for novel disease genes) and TR constraint scores
- Distribution of lengths/sd of lengths of repeats in CMH: repeats with low variability in the population can often come up as outliers if test sample is near the top of the distribution
- Consider also expanding to non-HPO genes and evaluate some of the top LPS ranked outliers to discover novel pathogenic repeat expansions in novel disease genes

> *Most* known pathogenic repeat expansions from our positive controls tend to be prioritized in the top 10 of ranked LPS outliers by these methods. They aren’t perfect but are a good primary analysis to identify large candidate expansions of pure reference motif repeats. You could also play around with sorting the report by proband allele length Z score to try to find non-reference repeat expansions, however, it’s very challenging to find pathogenic repeat expansions when considering allele length alone.

**Report columns**

| Column | Comment | Source | Example |
| --- | --- | --- | --- |
| CHROM | Chromosome # | TRGT | 1 |
| POS | Tandem repeat start position | TRGT | 1002390 |
| END | Tandem repeat end position | TRGT | 1002460 |
| TRID | Tandem repeat ID | TRGT | chr1_1002390_1002460_short |
| motif | Tandem repeat motif | TRGT | T |
| gene_name | Ensembl gene name | Ensembl, https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz | ISG15 |
| gene_id | Ensembl gene ID | Ensembl | ENSG00000187608 |
| gene_biotype | Gene type | Ensembl | protein_coding |
| feature | Gene feature | Ensembl | intron |
| omim_phenotype | OMIM gene phenotype associations | OMIM | Immunodeficiency 38, 616126 (3) |
| omim_inheritance | OMIM gene phenotype inheritance | OMIM | AR |
| HPO | HPO terms associated with proband phenotype | Phenotips | Seizures |
| gnomad_constraint_gene | Gene to which gnomAD LoF values are associated | GnomAD, https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv | ISG15 |
| lof.oe_ci.upper | LOEUF, upper bound of 90% confidence interval for o/e ratio for high confidence pLoF variants (lower values indicate more constrained, gnomAD recommends a threshold of \< 0.6 if needed) | gnomAD | 1.762 |
| lof.pLI | Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% o/e pLoF ratio; computed from high confidence pLoF gnomAD data) | gnomAD | 0.38914 |
| TR_constraint | Observed to expected ratio for length of the longest pure segment at TR locus ('TR constraint'). Pathogenic repeats exhibit higher LPS length variance than expected compared to rest of repeat catalog (TR_constraint \> 1, though protein-coding loci are under some constraint, OE \< 1). | Supplementary table 3 from 'Detailed tandem repeat allele profiling in 1,027 long-read genomes reveals genome-wide patterns of pathogenicity'. Note: the repeat loci don't line up exactly with the Adotto catalog, so I used an 80% reciprocal overlap to match the loci. Some loci therefore may not have a reported TR_constraint value. | 6.404 |
| ENCODE_promoter_ID | ENCODE cis-regulatory element promoter ID. Enter this ID into the search bar of the ENCODE SCREEN portal (https://screen.encodeproject.org/) for more details. | https://downloads.wenglab.org/cCREs/GRCh38-PLS.bed | EH38E1402406 |
| ENCODE_promoter_coord | Coordinates for promoters above. | Same as above | 1:181023144-181023324 |
| CMH_allele_len_range | Range of allele lengths in CMH controls | Generated value | (62, 77) |
| CMH_outlier_cutoff | Cutoff calculated from control allele lengths - mean plus standard deviation of the 95th quantile for TR allele length (resampled 100 times) | Generated value | 76 |
| CMH_allele_len_std | CMH allele length standard deviation. | Generated value | 254.23 |
| C4R_outlier_count | Number of expansions observed at a tandem repeat locus in the C4R cohort (139 individuals) relative to the CMH cohort (1436 individuals). | Generated value | 4 |
| C4R_outlier_samples | C4R sample IDs for samples in whom an expansion outlier was observed | Generated value | \<family_participant\> |
| C4R_outlier_allele_lens | C4R outlier allele lengths | Generated value | 1151;1220;1288;1775 |
| max_sample_allele_len | Maximum allele length observed across family | Generated value | 1192 |
| max_z_score_len | Maximum Z score of sample tandem repeat allele lengths in a family | Generated value | 3.17 |
| min_z_score_len_rank | Minimum Z score of sample tandem repeat allele lengths in a family relative to CMH, where 1 means it's the longest allele | Generated value | 38 |
| max_LPS | Maximum length of longest pure segment for this allele across family members (the longest uninterrupted repeat motif at a TR locus) | Generated value | 739 |
| min_LPS_rank | Minimum LPS of sample tandem repeat alleles in a family relative to CMH, where 1 means it has the longest LPS | Generated value | 2 |
| \<sample\>_LPS | Sample longest pure segment | TRGT-lps | 739 |
| num_outlier_samples | Number of samples with expanded allele (i.e. allele length \> cutoff and allele length Z score \>=3) | Generated value | 1 |
| \<sample\>_allele_len | Sample tandem repeat allele length. If the value is '.', TRGT failed to genotype this loci in this particular sample. | TRGT | 74 |
| \<sample\>_z_score_len | Sample tandem repeat allele length Z score (compared to CMH distribution). If the value is '.', TRGT failed to genotype this loci in this particular sample. | Generated value | 1.697 |
| \<sample\>_z_score_len_rank | Sample tandem repeat allele length Z score rank compared to CMH distribution. If the value is '.', TRGT failed to genotype this loci in this particular sample | Generated value | 38 |
| \<sample\>_AM | Sample tandem repeat average methylation | TRGT | 0.52 |
| \<sample\>_z_score_AM | Sample tandem repeat average methylation Z score (compared to control distribution). If missing, TRGT failed to assign methylation in this particular sample. | Generated value | 1.542 |
| \<sample\>_MP | Sample tandem repeat motif purity | TRGT | 0.511364 |
| \<sample\>_z_score_MP | Sample tandem repeat motif purity Z score (compared to control distribution) | Generated value | 3.174 |