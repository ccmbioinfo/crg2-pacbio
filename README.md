# Clinical research pipeline for exploring variants in PacBio HiFi whole genome (WGS) data (Hg38)

crg2-pacbio is a research pipeline aimed at discovering clinically relevant variants from PacBio HiFi whole genome sequence data in a family-based manner. crg2-pacbio uses Snakemake and Conda to manage jobs and software dependencies.

crg2-pacbio takes as input the following VCFs output by [PacBio's WGS pipeline](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL):
- [DeepVariant](https://github.com/google/deepvariant) joint-genotyped small variant VCF
- [pbsv](https://github.com/PacificBiosciences/pbsv) joint-genotyped structural variant VCF
- [TRGT](https://github.com/PacificBiosciences/trgt) VCFs, one per family member
- aligned BAM files for each sample

How to run the pipeline
1. Make a folder in a directory with sufficient space. Copy over the template files crg2-pacbio/samples.tsv, crg2-pacbio/units.tsv, crg2-pacbio/config.yaml, crg2-pacbio/crg2-pacbio.sh, crg2-pacbio/slurm_profile/slurm-config.yaml.  Note that 'slurm-config.yaml' is for submitting each rule as cluster jobs, so ignore this if not running on cluster.
```
mkdir NA12878
cp crg2-pacbio/samples.tsv crg2-pacbio/units.tsv crg2-pacbio/config.yaml crg2-pacbio/crg2-pacbio.sh crg2-pacbio/slurm_profile/slurm-config.yaml NA12878
cd NA12878
```
2. Set up pipeline run
- Reconfigure 'samples.tsv' and 'units.tsv' to reflect sample names and input files. Note that because several of the inputs are joint-genotyped, one row in the units.tsv file corresponds to one family, not one sample. 'units.tsv' must be configured with the path to the joint-genotyped deepvariant `small_variant_vcf`, the path to a directory containing per-sample genome-wide TRGT VCFs `trgt_vcf_dir`, and the joint-genotyped pbsv `pbsv_vcf`. The `trgt_vcf_dir` must ALSO contain the TRGT spanning BAMs to successfully run the denovo tandem repeat report. The 'samples.tsv' file should contain one row per family member, with the `sample` column corresponding to the family member ID and the `BAM` column corresponding to the path to the BAM file for that sample. 
units.tsv example:
```
family	platform	small_variant_vcf	trgt_vcf_dir	pbsv_vcf	trgt_pathogenic_vcf_dir
1042	PACBIO	/hpf/largeprojects/ccmbio/ccmmarvin_tcag_shared/PacBio/MEA26290/1042_TR0246.joint/1042_TR0246.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz	trgt	/hpf/largeprojects/ccmbio/ccmmarvin_tcag_shared/PacBio/MEA26290/1042_TR0246.joint/1042_TR0246.joint.GRCh38.pbsv.phased.vcf.gz
```
samples.tsv example:
```
sample	BAM
01	/hpf/largeprojects/ccmbio/ccmmarvin_tcag_shared/PacBio/MEA26290/1042_TR0246/01.m84090_240206_172710_s2.hifi_reads.bc2010.KL.GRCh38.aligned.haplotagged.bam
02	/hpf/largeprojects/ccmbio/ccmmarvin_tcag_shared/PacBio/MEA26290/1042_TR0247/02.m84090_240206_192642_s3.hifi_reads.bc2011.KL.GRCh38.aligned.haplotagged.bam
```
- Add paths to the HPO term file and pedigree file to config.yaml. 
- Do a dry run: add a `-n` flag to the Snakemake command in crg2-pacbio.sh. This will print out the rules that will be run, but not actually run them.
- If the dry run is successful, run the pipeline: `sbatch crg2-pacbio.sh`.
- If you just want to generate a single report, you can specify the report name in the Snakemake command in crg2-pacbio.sh. For example, to generate the repeat outlier report, add ` repeat_outliers/{family}.repeat.outliers.annotated.csv` to the Snakemake command.

Example:
```
#!/bin/bash
#SBATCH --job-name=crg2-pacbio
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out

SF=~/crg2-pacbio/Snakefile
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/snakemake"
SLURM=~/crg2-pacbio/slurm_profile/
CONFIG="config.yaml"

source /hpf/largeprojects/ccm_dccforge/dccdipg/Common/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --configfile ${CONFIG} --profile ${SLURM} repeat_outliers/{family}.repeat.outliers.annotated.csv

```

Small variant annotation and report: small_variants/coding/{family}/{family}.wes.regular.{date}.csv
- annotate variants with [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) and [vcfanno](https://github.com/brentp/vcfanno)
- generate a [gemini](https://github.com/arq5x/gemini) variant database
- generate a rare (less than 1% maximum gnomAD population AF) variant report with medium to high impact variants. Report generation scripts are derived from [cre](https://github.com/ccmbioinfo/cre), but copied to this repo for convenience. Report details are described [here](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/_layouts/15/Doc.aspx?sourcedoc=%7B7E4D3F4D-D83F-474C-A9CB-F59C2EB05C8A%7D&file=SNV_indel_report_November_2021.docx&action=default&mobileredirect=true); additional columns included are allele frequencies and counts from the [CoLoRSDB cohort](https://zenodo.org/records/11511513).

Small variant panel report: small_variants/panel/{family}/{family}.wgs.regular.{date}.csv, small_variants/panel-flank/{family}/{family}.wgs.regular.{date}.csv
- The annotation pipeline is the same as above, but only variants in a gene panel are considered. The panel report includes variants of any impact (i.e. it includes non-coding variants and intronic variants).
- To generate a gene panel from an HPO text file exported from PhenomeCentral or G4RD, add the HPO filepath to `config["run"]["hpo"]`.
- The first time you run the pipeline, you will also need to generate Ensembl and RefSeq gene files as well as an HGNC gene mapping file.
- Download and unzip Ensembl gtf: ```wget -qO- https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz  | gunzip -c > Homo_sapiens.GRCh38.112.gtf```
- Download and unzip RefSeq gff: ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz | gunzip -c > GRCh38_latest_genomic.gff```
- Download RefSeq chromosome mapping file: ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt```
- Run script to parse the above files: ```python scripts/clean_gtf.py --ensembl_gtf /path/to/Homo_sapiens.GRCh38.112.gtf --refseq_gff3 /path/to/GRCh38_latest_genomic.gff --refseq_assembly /path/to/GRCh38_latest_assembly_report.txt```
- Add the paths to the output files, Homo_sapiens.GRCh38.112.gtf_subset.csv and GRCh38_latest_genomic.gff_subset.csv, to the `config["gene"]["ensembl"]` and `config["gene"]["refseq"]` fields.
- You will also need the HGNC alias file: download this from https://www.genenames.org/download/custom/ using the default fields. Add the path this file to `config["gene"]["hgnc"]`.

Structural variant annotation and report: sv/{family}.pbsv.{date}.csv
- annotate variants using [SnpEFF](https://github.com/pcingola/SnpEff) and [AnnotSV](https://github.com/lgmgeo/AnnotSV)
- filter out variants under 50bp
- annotate SVs with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_SVs.py)
- a description of the report and fields can be found [here](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/_layouts/15/Doc.aspx?sourcedoc=%7B531618A5-B617-4444-B496-5A9D239C4B91%7D&file=PacBio_SV_report.docx&action=default&mobileredirect=true)

Repeat expansion outlier report: repeat_outliers/{family}.repeat.outliers.annotated.csv
- TRGT must have previously been run on each sample against the [937,122](https://zenodo.org/record/7987365#.ZHY9TOzMJAc) repeats originally released by the Genome in a Bottle tandem repeat benchmarking project
- this module is derived from the [find-outlier-expansions workflow](https://github.com/tandem-repeat-workflows/find-outlier-expansions/blob/main/find-outlier-expansions.ipynb) developed by Egor Dolzhenko (PacBio), Adam English (Baylor College of Medicine), Tom Mokveld (PacBio), Giulia Del Gobbo (CHEO), and Madeline Couse (SickKids)
- construct a repeat database from 98 HPRC samples and the individuals from the family of interest
- identify repeats with outlying size in family members
- annnotate repeat outliers with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_repeat_outliers.py)

De novo tandem repeat variant report: TRGT_denovo/{family}_{child}.TRGT.denovo.annotated.csv
- identify de novo tandem repeats in probands using [TRGT-denovo](https://github.com/PacificBiosciences/trgt-denovo) v0.2.0
- note that the `trgt_vcf_dir` must contain the TRGT spanning BAMs to successfully run the denovo tandem repeat report
- filter annotate de novo repeats with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_denovo_repeats.py)

Pathogenic repeat loci report: pathogenic_repeats/{family}.known.path.str.loci.csv
- genotype repeats per-sample using TRGTv1.0.0 against the [pathogenic repeat loci BED file](https://github.com/PacificBiosciences/trgt/blob/main/repeats/pathogenic_repeats.hg38.bed) provided by TRGT (the GIAB 937,122 catalog only contains 50/56 loci). Note that TRGT requires BAMs; these must be added to the samples.tsv file
- merge sample VCFs into multi-sample family VCF
- annotate repeat loci with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_path_str_loci.py)


