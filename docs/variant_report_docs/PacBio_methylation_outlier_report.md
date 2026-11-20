# Methylation outlier report

Madeline Couse, Giulia del Gobbo, Christian Marshall, Kristin Kernohan

**Version 2025-09**

## Pipeline

Each case is run a background of controls (for example, all C4R long-read genomes or all TG long-read genomes).

The methylation workflow is as follows:

1.  Extract CpG methylation probabilities from PacBio HiFi BAMs using [pb-cpg-tools](https://github.com/PacificBiosciences/pb-CpG-tools) v2.3.2 with a coverage threshold of 4.

2.  Generate methylation profiles for 200bp tiles across the genome for each sample using [Methbat](https://github.com/PacificBiosciences/MethBat/tree/main) v0.13.3.

3.  Calculate average coverage at each 200bp tile using mosdepth v0.3.11.

4.  Build a MethBat background profile using the control samples.

5.  Identify methylation outliers in sample of interest versus controls using MethBat.

6.  Filter, merge, and annotate methylation outliers with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_methylation_outliers.py)

    1.  Filter out regions where ‘compare_label’ column (described in table below) is ‘Uncategorized’ (the region does not match any of HyperASM, HypoASM, HypoMethylated, or HyperMethylated labels) or ‘InsufficientData’ (the dataset does not have CpGs or the number of datasets in the background is not sufficient).

    2.  Merge adjacent methylation outlier tiles, counting the number of contributing tiles in the merged tile (‘tile_count’ column). Aggregation functions applied to per-tile MethBat annotations while merging adjacent tiles are detailed in Table 1.

    3.  Annotate methylation outliers against nearby variants

        1.  SNVs and indels within 50bp of outlier tile: less than 1% gnomAD population maximum allele frequency

        2.  SVs and CNVs within 1000 bp of outlier tile. No frequency filter as it is difficult to choose an appropriate filter.

        3.  Tandem repeats with Z-score greater than or equal to 3 within 1000 bp of outlier tile

    4.  Annotate methylation outliers against genes, OMIM phenotypes, HPO terms, constraint scores, and regulatory regions.

| **Field** | **Aggregation Method** |
| --- | --- |
| POS | minimum |
| END | maximum |
| summary_label | joined with ';' as separator |
| comparison_label | joined with ';' as separator |
| category_pop_count | maximum |
| category_pop_freq | maximum |
| asm_fishers_pvalue | minimum |
| mean_hap1_methyl | mean |
| mean_hap2_methyl | mean |
| mean_meth_delta | maximum |
| mean_abs_meth_delta_zscore | minimum if minimum zscore \< 0, maximum otherwise |
| mean_combined_methyl | mean |
| mean_combined_methyl_zscore | minimum if minimum zscore \< 0, maximum otherwise |
| num_phased_cpgs | summed |
| num_unphased_cpgs | summed |
| num_partial_cpgs | summed |

Table 1. Aggregation methods applied to MethBat annotations when merging adjacent outlier tiles.

## Recommended analysis strategy 

The methylation outliers are already sorted in descending order by the absolute Z-score values. It’s not a bad idea to filter leniently (e.g. only filter on population count, presence in a gene) and do a first pass on those variants to see if there are any rare outliers with large Z-scores. If nothing stands out, filter further. See Table 2 for filtering recommendations. We tested this pipeline on six positive controls from Children’s Mercy Hospital and found that after filtering on population count (\<2), gene, and maximum number of CpGs (\>=10), the validated methylation outliers were ranked first when sorted by Z-score (Z-scores ranged from 13-34.5).

Tips for interpretation:

- ‘tile_count’: - multiple adjacent tiles with the same methylation outlier call is compelling for a large-impact change. In our positive controls from CMH, top outliers related to diagnostic typically were associated with outlier calls in multiple adjacent tiles
- ‘nearby_variant’: most DNA methylation changes caused by an underlying genetic variant are within a short distance of the given variant. This column is helpful to highlight whether any rare small variants, outlier TRs, or any (\*note not filtered for rarity) SVs are nearby to the outlier tile(s). Use this to help prioritize outlier tiles, but consider also that more long-range epigenetic effects are possible, or that in very rare cases the mechanism of disease may be epigenetic rather than genetic
- ‘GREENDB’ columns: help interpret whether methylation outliers overlap annotated regulatory elements
- In general, use gene-level annotations as you would for interpretation of other variant types (e.g. consider OMIM phenotypes, gene constraint scores, etc.)

The selection contains two HTML tables. Both are converted below, with the caption and footnote between them left as they are.

| Field | Filter |
| --- | --- |
| tile_count | \> 1 |
| category_pop_count | ≤ 3. Can modify depending on expected inheritance and number of affected per family |
| gene_id | != "." |
| max_num_cpgs* | Play around with this. Start with no filter to see if there are any top hits of interest, and go as high as \>=10 to reduce the number of candidates considerably. A variant can result in a methylation change at a single CpG site in a region, but I'm not sure that this would correlate with the sort of large effect sizes we expect to see with rare variants. |
| nearby_variant | != ".' In theory, we're looking for methylation changes associated with rare variation, so we'd expect to see a rare variant nearby (and can make a more convincing argument that the methylation change is rare if it's linked to a very rare variant). That said, there can certainly be methylation changes that are the result of changes elsewhere in the genome that are of interest. And, although not as common, a rare disease may even have an epigenetic basis as the result of some embryonic insult. |

Table 2. Filtering recommendations

\*The max_num_cpgs column refers to the number of CpG sites in an outlier region that MethBat used to compute methylation statistics. This number may not represent the actual number of CpGs in a region, because pb-cpg-tools requires coverage of at least five reads to produce a methylation probability for a CpG site. So, if coverage is low in a region, the methylation probability is not called and therefore that site is not included in MethBat's calculations.

**Report columns**

| Column | Comment | Source | Example |
| --- | --- | --- | --- |
| CHROM | Chromosome # | Calculated | 1 |
| POS | Tile start position | Calculated | 1002390 |
| END | Tile end position | Calculated | 1002460 |
| tile_count | Number of adjacent outlier tiles | Calculated | 2 |
| summary_label | **Methylated** - indicates that the combined CpG methylation had a high average methylation rate (by default, \>=80%). **Unmethylated** - indicates that the combined CpG methylation had a low average methylation rate (by default, \<=20%). **AlleleSpecificMethylation** - indicates that a sufficient fraction of the CpGs were phased (by default, \>=50%) and that ASM was detected through both a significant Fisher's exact test (by default, p \<= 0.01) and difference in mean methylation for haplotypes 1 and 2 (by default, \>= 40% methylation delta) | MethBat | AlleleSpecificMethylation |
| compare_label | High level summary label of the comparison of this dataset to the background. **HypoASM** - if the region shows significantly less ASM in the dataset than the background; this option has lower priority than the following options. **HyperASM** - if the region shows significantly more ASM in the dataset than the background; the summary_label must be "AlleleSpecificMethylation". **HypoMethylated** - if the region shows significantly less methylation in the dataset than the background; the summary_label must be "Unmethylated". **HyperMethylated** - if the region shows significantly more methylation in the dataset than the background; the summary_label must be "Methylated" | MethBat | HyperASM |
| gene_name | Ensembl gene name | Ensembl, https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz | ISG15 |
| gene_id | Ensembl gene ID | Ensembl | ENSG00000187608 |
| gene_biotype | Gene type | Ensembl | protein_coding |
| feature | Gene feature | Ensembl | intron |
| nearby_variant | Closest rare small variant (within 50bp), or SV, or CNV, or TR (Z-score \> 3 relative to population mean in the Children's Mercy Hospital GA4K cohort) within 1000bp. 1. Rare (less than 0.01 gnomAD_af_popmax) SNVs and indels: small_var:{chr}:{pos}-{end}-{ref}-{alt}-{genotype}-{gnomad_AF_popmax}-{distance_from_tile}. 2. SVs: SV: {chr}:{pos}-{end}-{svtype}-{genotype}-{cmh_max_AF or CoLoRSDB_max_AF}-{distance_from_tile}. 3. CNVs: CNV: {chr}:{pos}-{end}-{svtype}-{genotype}-{distance_from_tile}. 4. Tandem repeats with allele length Z-score \>=3: TR:chr}:{pos}-{end}-{motif}-{allele length Z-score}-{distance_from_tile} | calculated | TR:15:25389858-25389952-ACAG-TR-zscore:6.323-dist:459 |
| closest_exon_gene | Gene name (or ID, if name is missing) associated with closest exon | Ensembl | ISG15 |
| omim_phenotype | OMIM gene phenotype associations | OMIM | Immunodeficiency 38, 616126 (3) |
| omim_inheritance | OMIM gene phenotype inheritance | OMIM | AR |
| HPO | HPO terms associated with proband phenotype | Phenotips | Seizures |
| gnomad_constraint_gene | Gene to which gnomAD LoF values are associated | GnomAD, | ISG15 |
| lof.oe_ci.upper | LOEUF, upper bound of 90% confidence interval for o/e ratio for high confidence pLoF variants (lower values indicate more constrained, gnomAD recommends a threshold of \< 0.6 if needed) | gnomAD | 1.762 |
| lof.pLI | probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% o/e pLoF ratio; computed from high confidence pLoF gnomAD data) | gnomAD | 0.38914 |
| category_pop_count | The number of datasets in the background population with the same summary_label | MethBat | 2 |
| category_pop_freq | The percentage of datasets (excluding NoData) in the background population with the same summary_label | MethBat | 0.001 |
| asm_fishers_pvalue | The raw p-value from a Fisher's exact test comparing the two haplotypes and the number of reads that are methylated/unmethylated | MethBat | 0 |
| mean_hap1_methyl | The mean (average) methylation ratio for CpGs on an individual haplotype 1 | MethBat | 0.53 |
| mean_hap2_methyl | The mean (average) methylation ratio for CpGs on an individual haplotype 2 | MethBat | 0.02 |
| mean_meth_delta | The difference in mean methylation ratios between the two haplotypes, mean_hap2_methyl - mean_hap1_methyl | MethBat | -0.51 |
| mean_abs_meth_delta_zscore | If a background region file was provided, this value is the Z-score comparing abs(mean_meth_delta) against the background profile; positive values indicate more evidence of ASM in this dataset relative to the population, negative values indicate less | MethBat | 10.54 |
| mean_combined_methyl | The mean (average) combined methylation ratio; "combined" here indicates that phasing (i.e. haplotypes) is not considered | MethBat | 0.04 |
| mean_combined_methyl_zscore | If a background region file was provided, this value is the Z-score comparing mean_combined_methyl against the background profile; positive values indicate this dataset is hyper-methylated relative to the the population, negative values indicate hypo-methylation | MethBat | 0.51 |
| max_abs_meth_delta_zscore | Maximum of absolute mean_abs_meth_delta_zscore and mean_combined_methyl_zscore Z scores | Calculated | 10.54 |
| num_phased_cpgs | The number of CpGs in the region with haplotagged reads on only one haplotype | MethBat | 12 |
| num_partial_cpgs | The number of CpGs in the region with haplotagged reads on only one haplotype | MethBat | 0 |
| num_unphased_cpgs | The number of CpGs in the region with no haplotagged reads | MethBat | 0 |
| max_num_cpgs | Maximum of num_phased_cpgs and num_unphased_cpgs | Calculated | 12 |
| coverage | Average coverage in tile, calculated by mosdepth | mosdepth | 40.06 |
| GREENDB_reg_region | Regulatory element annotations from GREENDB ( are using the GREEN-DB bed file) | GREENDB | Enhancer;promoter |
| GREENDB_source | Regulatory region annotation source | GREENDB | EnsemblRegBuild,ENCODE-HMM;;SegWey |
| GREENDB_closest gene | Gene closest to regulatory element | GREENDB | WBP4;WBP4 |
| GREENDB_controlled_genes | Gene(s) predicted to be controlled by regulatory element | GREENDB | WBP4;WBP4 |