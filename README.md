# Clinical research pipeline for exploring variants in PacBio HiFi whole genome (WGS) data (Hg38)

crg2-pacbio is a research pipeline aimed at discovering clinically relevant variants from PacBio HiFi whole genome sequence data in a family-based manner. crg2-pacbio uses Snakemake and Conda to manage jobs and software dependencies.

crg2-pacbio takes as input the following VCFs output by [PacBio's WGS pipeline](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL):
- [DeepVariant](https://github.com/google/deepvariant) joint-genotyped small variant VCF
- [pbsv](https://github.com/PacificBiosciences/pbsv) joint-genotyped structural variant VCF
- [TRGT](https://github.com/PacificBiosciences/trgt) VCFs, one per family member

Small variant annotation and report
- annotate variants with [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) and [vcfanno](https://github.com/brentp/vcfanno)
- generate a [gemini](https://github.com/arq5x/gemini) variant database
- generate a rare (less than 1% maximum gnomAD population AF) variant report with medium to high impact variants. Report generation scripts are derived from [cre](https://github.com/ccmbioinfo/cre), but copied to this repo for convenience. Report details are described [here](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/_layouts/15/Doc.aspx?sourcedoc=%7B7E4D3F4D-D83F-474C-A9CB-F59C2EB05C8A%7D&file=SNV_indel_report_November_2021.docx&action=default&mobileredirect=true); additional columns included are allele frequencies and counts from the [CoLoRSDB cohort](https://zenodo.org/records/11511513).

Small variant panel report
- The annotation pipeline is the same as above, but only variants in a gene panel are considered. The panel report includes variants of any impact (i.e. it includes non-coding variants and intronic variants).
- To generate a gene panel from an HPO text file exported from PhenomeCentral or G4RD, add the HPO filepath to `config["run"]["hpo"]`.
- The first time you run the pipeline, you will also need to generate Ensembl and RefSeq gene files as well as an HGNC gene mapping file.
- Download and unzip Ensembl gtf: ```wget -qO- https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz  | gunzip -c > Homo_sapiens.GRCh38.112.gtf```
- Download and unzip RefSeq gff: ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz | gunzip -c > GRCh38_latest_genomic.gff```
- Download RefSeq chromosome mapping file: ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt```
- Run script to parse the above files: ```python scripts/clean_gtf.py --ensembl_gtf /path/to/Homo_sapiens.GRCh38.112.gtf --refseq_gff3 /path/to/GRCh38_latest_genomic.gff --refseq_assembly /path/to/GRCh38_latest_assembly_report.txt```
- Add the paths to the output files, Homo_sapiens.GRCh38.112.gtf_subset.csv and GRCh38_latest_genomic.gff_subset.csv, to the `config["gene"]["ensembl"]` and `config["gene"]["refseq"]` fields.
- You will also need the HGNC alias file: download this from https://www.genenames.org/download/custom/ using the default fields. Add the path this file to `config["gene"]["hgnc"]`.

Structural variant annotation and report
- annotate variants using [SnpEFF](https://github.com/pcingola/SnpEff) and [AnnotSV](https://github.com/lgmgeo/AnnotSV)
- filter out variants under 50bp
- annotate SVs with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_SVs.py)
- a description of the report and fields can be found [here](https://sickkidsca.sharepoint.com/:w:/r/sites/thecenterforcomputationalmedicineworkspace/_layouts/15/Doc.aspx?sourcedoc=%7B531618A5-B617-4444-B496-5A9D239C4B91%7D&file=PacBio_SV_report.docx&action=default&mobileredirect=true)

Repeat expansion outlier report
- TRGT must have previously been run on each sample against the [937,122](https://zenodo.org/record/7987365#.ZHY9TOzMJAc) repeats originally released by the Genome in a Bottle tandem repeat benchmarking project
- this module is derived from the [find-outlier-expansions workflow](https://github.com/tandem-repeat-workflows/find-outlier-expansions/blob/main/find-outlier-expansions.ipynb) developed by Egor Dolzhenko (PacBio), Adam English (Baylor College of Medicine), Tom Mokveld (PacBio), Giulia Del Gobbo (CHEO), and Madeline Couse (SickKids)
- construct a repeat database from 98 HPRC samples and the individuals from the family of interest
- identify repeats with outlying size in family members
- annnotate repeat outliers with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_repeat_outliers.py)

Pathogenic repeat loci report
- genotype repeats per-sample using TRGTv1.0.0 against the [pathogenic repeat loci BED file](https://github.com/PacificBiosciences/trgt/blob/main/repeats/pathogenic_repeats.hg38.bed) provided by TRGT (the GIAB 937,122 catalog only contains 50/56 loci). Note that TRGT requires BAMs; these must be added to the samples.tsv file
- merge sample VCFs into multi-sample family VCF
- annotate repeat loci with [custom Python script](https://github.com/ccmbioinfo/crg2-pacbio/blob/master/scripts/annotate_path_str_loci.py)


