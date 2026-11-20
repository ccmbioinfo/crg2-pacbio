# PacBio TRGT genotypes for known pathogenic tandem repeats

Madeline Couse, Christian Marshall, Kristin Kernohan, Giulia Del Gobbo

Last updated: October 2025

[TRGT](https://github.com/PacificBiosciences/trgt) v2.1.0 is a tandem repeat genotyper that allows targeted genotyping of tandem repeats from PacBio HiFi data. We genotype a catalog of repeats derived from the union of the [STRchive TRGT definitions](https://strchive.org/_astro/STRchive-disease-loci.hg38.TRGT.B2FLLAlV.bed) and the [bed file provided by TRGT](http://null/#L2) file that together comprise the coordinates and motifs of 77 repeat loci. Expansions in the majority of these loci are associated with disease. We run TRGT using this catalog on PacBio genomes (individually, as joint-genotyping is not a feature), and then merge the VCFs by family. Lastly, we create a CSV report using the resulting TRGT VCF. We pulled disease thresholds from the [STRchive BED file for filtering genomic data](https://strchive.org/_astro/STRchive-disease-loci.hg38.CE2vK2zA.bed), as well as the literature. Descriptions for the report columns are listed in the table below.

Suggestions for filtering and interpretation

- Filter out loci that are homozygous reference in all family members: ALT != ‘homozygous_ref’
- Filter for loci that have a putative expansion, or complex loci: DISEASE_PREDICTION contains ‘True’ OR ‘None’

NB:

For sample motif_counts with an underscore (e.g. 12_0\|14_0), we do not compute a disease prediction. Instead, this column is set to ‘None’ for all samples. Whether repeats are beyond the pathogenic threshold will have to be manually interpreted by analysts based on the DISEASE_THRESHOLD and motif_count and/or allele_length columns.

At the *VWA1* locus, any deviation from two repeat copies is thought to be pathogenic, i.e. contractions or expansions are pathogenic (STRchive).

For detailed descriptions of repeat loci, including descriptions, pathogenic expansion repeat ranges, prevalence, age of onset, references and more, please refer to the [STRchive loci page](https://strchive.org/loci/).

## Column descriptions

| Column | Source | Example |
| --- | --- | --- |
| CHROM | TRGT VCF | chr10 |
| POS | TRGT VCF | 93702523 |
| END | TRGT VCF | 93702547 |
| REF | TRGT VCF | CCGC |
| ALT | TRGT VCF | CCACCGCCGCCGCCGCCGCCGCCGCCGC\|C\| *This may be 'homozygous_ref' if the locus is homozygous reference in all family members |
| GENE | TRGT VCF | FRA10AC1 |
| DISORDER | /hpf/largeprojects/ccm_dccforge/dccdipg/Common/annotation/ExpansionHunter/tandem_repeat_disease_loci_v1.1.tsv | FRA10A: No confirmed disease association |
| MOTIFS | TRGT VCF | CCG |
| STRUC | TRGT VCF | (CCG)n |
| DISEASE_THRESHOLD | /hpf/largeprojects/ccm_dccforge/dccdipg/Common/annotation/ExpansionHunter/tandem_repeat_disease_loci_v1.1.tsv | 200.0 |
| DISEASE_PREDICTION | Calculated | False\|False\|False\|False |
| \<SAMPLE_ID\>_motif_count | TRGT VCF | 8\|8 |
| \<SAMPLE_ID\>_allele_length | TRGT VCF | 29\|32 |
| \<SAMPLE_ID\>_allele_CI | TRGT VCF | 28-30\|31-34 |
| \<SAMPLE_ID\>_spanning_read_support | TRGT VCF | 16\|11 |
| \<SAMPLE_ID\>_methylation | TRGT VCF | 0.259\|0.839 |
