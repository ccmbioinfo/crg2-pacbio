# Clinical research pipeline for exploring variants in PacBio HiFi whole genome (WGS) data

Contributors (in alphabetical order)
- Madeline Couse, SickKids
- Giulia Del Gobbo, CHEO
- Egor Dolzhenko, PacBio (repeat expansion outliers)
- Anjali Jain, SickKids

crg2-pacbio is a research pipeline aimed at discovering clinically relevant variants from PacBio HiFi whole genome sequence data in a family-based manner. crg2-pacbio uses Snakemake and Conda to manage jobs and software dependencies.

crg2-pacbio takes as input the following VCFs output by [PacBio's WGS pipeline](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL):
- [DeepVariant](https://github.com/google/deepvariant) joint-genotyped small variant VCF
- [pbsv](https://github.com/PacificBiosciences/pbsv) joint-genotyped structural variant VCF
- [TRGT](https://github.com/PacificBiosciences/trgt) VCFs, one per family member

Small variant annotation and report
- annotate variants with [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) and [vcfanno](https://github.com/brentp/vcfanno)
- generate a [gemini](https://github.com/arq5x/gemini) variant database
- generate a rare (less than 1% maximum gnomAD population AF) variant report with medium to high impact variants. Report generation scripts are derived from [cre](https://github.com/ccmbioinfo/cre), but copied to this repo for convenience. Report details are described [here](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/_layouts/15/Doc.aspx?sourcedoc=%7B7E4D3F4D-D83F-474C-A9CB-F59C2EB05C8A%7D&file=SNV_indel_report_November_2021.docx&action=default&mobileredirect=true); additional columns included are allele frequencies and counts from the [Children's Mercy Hospital (CMH, 502 genomes) cohort](https://github.com/ChildrensMercyResearchInstitute/GA4K/tree/main/pacbio_sv_vcf) and the [Human Pangenome Reference Consortium](https://humanpangenome.org/) (HPRC, 40 genomes).

Structural variant annotation and report
- annotate variants using [SnpEFF](https://github.com/pcingola/SnpEff) and [AnnotSV](https://github.com/lgmgeo/AnnotSV)
- filter out variants under 50bp
- annotate SVs with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_SVs.py)
- a description of the report and fields can be found [here](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/_layouts/15/Doc.aspx?sourcedoc=%7B531618A5-B617-4444-B496-5A9D239C4B91%7D&file=PacBio_SV_report.docx&action=default&mobileredirect=true)

Repeat expansion outlier report
- TRGT must have previously been run on each sample against the [937,122](https://zenodo.org/record/7987365#.ZHY9TOzMJAc) repeats originally released by the Genome in a Bottle tandem repeat benchmarking project
- this module is derived from the [find-outlier-expansions workflow](https://github.com/tandem-repeat-workflows/find-outlier-expansions/blob/main/find-outlier-expansions.ipynb) developed by Egor Dolzhenko, Giulia Del Gobbo, and Madeline Couse 
- construct a repeat database from 98 HPRC samples and the individuals from the family of interest
- identify repeats with outlying size in family members
- annnotate repeat outliers with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_repeat_outliers.py)




