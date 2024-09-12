# tdb multi-sample TRGT merge and tandem repeat statistics scripts

This repository contains scripts for creating, merging, and analyzing [tdb](https://github.com/ACEnglish/tdb) (tandem repeat database) files.
Thank you to Adam English (Baylor College of Medicine) for his assistance with multi-level tdb merging. 

## Scripts

### 1. create_tdb.sh

This script generates TDB files for each sample from TRGT VCF files.

#### Usage
./create_tdb.sh <vcf_dir>
Where `<vcf_dir>` is the directory containing the TRGT VCF files.

#### Process

1. Creates a list of VCF files.
2. Generates job scripts for each sample.
3. Submits jobs to create tdb files for each sample.

### 2. merge_tdbs.sh

This script performs an intermediate merge. So for example if we have 100 tdbs from 100 TRGT VCFs, we may merge tdbs 10 at a time, resulting in 10 tdbs.

#### Usage
./merge_tdbs.sh <tdb_dir>
Where `<tdb_dir>` is the directory containing the tdb files to be merged.

#### Process

1. Creates necessary directories for jobs, logs, and output.
2. Generates a job script for the intermediate merges.
3. Submits the jobs to a SLURM scheduler.

### 3. bigmerge_tdbs_level2_to_final.sh

This script performs the final merge of tdb files for projects requiring two merge levels (e.g., TCAG). 
If you have many VCFs (e.g. >1000), you may wish to use bigmerge_tdbs_level3.sh.

#### Usage
./bigmerge_tdbs_level2_to_final.sh <tdb_dir>
Where `<tdb_dir>` is the directory containing the tdb files to be merged.

#### Process

1. Creates necessary directories for jobs, logs, and output.
2. Generates a job script for the final merge.
3. Submits the jobs to a SLURM scheduler.

### 4. bigmerge_tdbs_level3.sh

This script performs the final merge of tdb files for projects requiring three merge levels (e.g., CMH).

#### Usage
bash
./bigmerge_tdbs_level3.sh <tdb_dir>

Where `<tdb_dir>` is the directory containing the level 2 tdb files to be merged.

#### Process

1. Removes any existing final merge tdb.
2. Creates a directory for the final tdb.
3. Generates a job script for the final merge.
4. Submits the job to a SLURM scheduler.

### 5. tdb_to_allele_dist.sh

This script runs a Python script to convert tdb data to allele distribution data.

#### Usage

Submit this script as a SLURM job:
sbatch tdb_to_allele_dist.sh

#### Process

1. Loads the Anaconda module.
2. Activates the tdb environment.
3. Runs the `tdb_to_allele_dist.py` Python script. This script calculates tandem repeat allele length and methylation distributions from a tdb. The output has the following format:

| LocusID | trid | AL_range | AL_95_mean | AL_95_std | AM_mean | AM_std |
|---------|------|----------|------------|-----------|---------|--------|
| 0 | chr1_83766_84081 | [284, 401] | 350.0 | 25.5 | 0.657 | 0.0813 |

## Requirements

- SLURM job scheduler
- Conda environment with tdb tools installed
- Python 3 with necessary libraries for tdb analysis


## Notes

- The scripts use a Conda environment located at `/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/envs/tdb/`.
- The merge processes use varying amounts of memory and CPU threads depending on the level of merge.
- The `tdb_to_allele_dist.sh` script requires may require a large amount of RAM (e.g. 300GB) depending on the number of samples in the tdb.

