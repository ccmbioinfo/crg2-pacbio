# *De novo* tandem repeat mutations report

Madeline Couse, Giulia del Gobbo, Christian Marshall, Kristin Kernohan

**Version 2024-10**

The *de novo* tandem repeat (TR) expansion workflow uses [TRGT-denovo](https://github.com/PacificBiosciences/trgt-denovo/tree/main) ([link to preprint](https://doi.org/10.1101/2024.07.16.600745)) to identify *de novo* TR expansions, contractions, and compositional changes in parent-child trios. TRGT-denovo integrates spanning reads from BAMs with TRGT genotyping results to reduce the number of false positive *de novo* TR mutations. As with the tandem repeat outlier workflow, we run TRGT (v0.5.0) on the [adotto repeat catalog](https://zenodo.org/record/7987365#.ZHY9TOzMJAc) (937,122 repeat loci). Then, TRGT-denovo gathers evidence for possible de novo tandem repeat mutations from the TRGT VCF and BAM.

**Filters used to define candidate *de novo* TRs**

The TRGT-denovo (v0.2.0) raw output is filtered by the criteria suggested by the creators of TRGT-denovo in this [Jupyter Notebook](https://github.com/PacificBiosciences/trgt-denovo/blob/main/scripts/python/trio_analysis.ipynb).

| **Field** | **Meaning** | **Filter** |
| --- | --- | --- |
| denovo_coverage | Number of child reads supporting a de novo allele | \>= 5 |
| allele_ratio | Proportion of reads in an allele having *de novo* signal (*de novo* reads / all reads belonging to that allele), where this is ideally close to 1 to indicate that there is little ambiguity of there being a *de novo* call (lower ratios may indicate many possible causes, e.g., mosaicism, somatic instability, or stutter) | \>= 0.7 |
| child_ratio | Proportion of reads at a locus having *de novo* signal (*de novo* reads / all reads covering a locus), which ideally approaches 0.5 | 0.3 \<= child ratio \<= 0.7 unless variant is on the X chromosome |
| mother_dropout_prob | Allele dropout rate for reads coming from the mother | \< 0.005 |
| father_dropout_prob | Allele dropout rate for reads coming from the father | \< 0.005 |
| mother_dropout | Coverage cut-off dropout detection using HP tags from phasing tools in mother; possible values: Full dropout (FD), Haplotype dropout (HD), Not (N) | N |
| father_dropout | Coverage cut-off dropout detection using HP tags from phasing tools in father; possible values: Full dropout (FD), Haplotype dropout (HD), Not (N) | N |
| mean_diff_mother \> 7 | Score difference between de novo and maternal reads; lower values indicate greater similarity | \>7 |
| mean_diff_father \>7 | Score difference between de novo and paternal reads; lower values indicate greater similarity | \>7 |

**Control tandem repeat allele length distributions**

We calculated the TR allele length distributions for the [adotto repeat catalog](https://zenodo.org/record/7987365#.ZHY9TOzMJAc) (937,122 repeat loci) from 1,436 individuals from the Children’s Mercy GA4K cohort. We calculated the allele length distributions for the short allele and long alleles separately, i.e. we aggregated the shortest allele for all samples and calculated sample statistics and did the same for the longest allele for all samples. This may be useful in autosomal recessive TR expansion inheritance; a TR expansion may not be an outlier if the allele is present in the population, but if we consider the short and long alleles separately, the 'short' allele should be picked up as an outlier if there are few or no homozygous carriers in the population. The cutoffs for the short and long alleles, used to determine repeat alleles expanded relative to the control population, are calculated as the mean plus standard deviation of the 95th quantile for TR allele length (resampled 100 times).

**Recommended filter for analysis**

C4R_outlier_count \< 3

outlier != False

gene_name != ‘.’

denovo_status != Y:- or Y:= (allele is not a contraction or substitution)

**Report columns**

| Column | Comment | Source | Example |
| --- | --- | --- | --- |
| CHROM | Chromosome # | TRGT | 11 |
| POS | Tandem repeat start position | TRGT | 47633475 |
| END | Tandem repeat end position | TRGT | 47633552 |
| TRID | Tandem repeat ID | TRGT | chr11_47633475_47633552_T_long_allele |
| gene_name | Ensembl gene name | Ensembl, https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz | ISG15 |
| gene_id | Ensembl gene ID | Ensembl | ENSG00000187608 |
| gene_biotype | Gene type | Ensembl | protein_coding |
| segdup | Segmental duplication | PacBio resources | chr1:1,464,080-1,496,342 |
| omim_phenotype | OMIM gene phenotype associations | OMIM | Immunodeficiency 38, 616126 (3) |
| omim_inheritance | OMIM gene phenotype inheritance | OMIM | AR |
| HPO | HPO terms associated with proband phenotype | Phenotips | Seizures |
| gnomad_constraint_gene | Gene to which gnomAD LoF values are associated | GnomAD, https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv | ISG15 |
| lof.oe_ci.upper | LOEUF, upper bound of 90% confidence interval for o/e ratio for high confidence pLoF variants (lower values indicate more constrained, gnomAD recommends a threshold of \< 0.6 if needed) | gnomAD | 0.689 |
| lof.pLI | probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% o/e pLoF ratio; computed from high confidence pLoF gnomAD data) | gnomAD | 0.14085 |
| feature | Gene feature | Ensembl | intron |
| child_AL | Child allele lengths at this locus | TRGT | 67,117 |
| mother_AL | Mother allele lengths at this locus | TRGT | 67,102 |
| father_AL | Father allele lengths at this locus | TRGT | 67,105 |
| outlier | Child *de novo* TR allele length is greater than or equal to the control range estimated from the CMH cohort | Generated value | TRUE |
| z_score_len | Child *de novo* TR allele length z score (compared to CMH control distribution). | Generated value | 1.864 |
| control_range | Range of allele lengths in controls | Generated value | (66, 132) |
| cutoff | Cutoff calculated from control allele lengths - mean plus standard deviation of the 95th quantile for TR allele length (resampled 100 times) | Generated value | 111.94 |
| C4R_outlier_count | Number of expansions observed at a tandem repeat locus in the C4R cohort (139 individuals) relative to the Children’s Mercy cohort (1436 individuals). | Generated value | 2 |
| C4R_outlier_samples | C4R sample IDs for samples in whom an expansion outlier was observed | Generated value | 1074_CH0164 |
| denovo_coverage | Number of child reads supporting a de novo allele compared to parental data | TRGT-denovo | 8 |
| allele_coverage | Number of child reads mapped to this specific allele | TRGT-denovo | 8 |
| allele_ratio | Proportion of reads in an allele having *de novo* signal (*de novo* reads / all reads belonging to that allele), where this is ideally close to 1 to indicate that there is little ambiguity of there being a *de novo* call (lower ratios may indicate many possible causes, e.g., mosaicism, somatic instability, or stutter) | TRGT-denovo | 1.0 |
| child_coverage | Total number of child reads at this site | TRGT-denovo | 18 |
| child_ratio | Proportion of reads at a locus having *de novo* signal (*de novo* reads / all reads covering a locus), which ideally approaches 0.5 | TRGT-denovo | 0.444 |
| mean_diff_father | Score difference between de novo and paternal reads; lower values indicate greater similarity | TRGT-denovo | 13.75 |
| mean_diff_mother | Score difference between de novo and maternal reads; lower values indicate greater similarity | TRGT-denovo | 33.75 |
| allele_origin | Inferred origin of the allele based on alignment; possible values: {F:{1,2,?}, M:{1,2,?}, ?}. F and M denote father and mother respectively. The associated {1, 2, ?} values denote the first or second allele from either parent or ? when this cannot be derived unambiguously. Lastly a ? denotes an allele for which parental origin cannot be determined unambiguously | TRGT-denovo | F:2 |
| denovo_status | ‘Y’ indicates a denovo allele is present (all ‘X’ are filtered out). If parental origin can be determined without ambiguity the allele sequences can be compared directly such that the de novo type can be established as + (expansion), - (contraction), or = (substitution) | TRGT-denovo | Y:+ |
| per_allele_reads_father | Number of reads partitioned per allele in the father (allele1, allele2) | TRGT-denovo | 11,11 |
| per_allele_reads_mother | Number of reads partitioned per allele in the mother (allele1, allele2) | TRGT-denovo | 19,13 |
| per_allele_reads_child | Number of reads partitioned per allele in the child (allele1, allele2) | TRGT-denovo | 10,8 |
| father_MC | TRGT VCF motif counts for this locus in the father | TRGT-denovo | 67,105 |
| mother_MC | TRGT VCF motif counts for this locus in the mother | TRGT-denovo | 67,102 |
| child_MC | TRGT VCF motif counts for this locus in the child | TRGT-denovo | 67,117 |
| min_mean_diff | Minimum of mean_diff_father and mean_diff_mother; approximates the relative size of the *de novo* event | TRGT-denovo | 13.75 |
