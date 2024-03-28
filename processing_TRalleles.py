#Developing a Database that tells us the min&max repeat sizes in samples for each locus
#Columns: repeat	min_size	max_size	min_size_sample	max_size_sample (5 columns)

from pathlib import Path
import gzip
import itertools
import numpy as np
import pandas as pd
from collections import namedtuple
import argparse 
from datetime import date


def get_repeat_recs(path):
    '''
    Parse through file and extract the repeat information (sample, alleles' size) along with their chrom location
    '''
    RepeatRec = namedtuple("RepeatRec", "sample short_allele long_allele")
    def parse_alleles(group):
        alleles = list(line.decode("utf8").split() for line in group)
        alleles = [(rec[1], rec[2].split(",")) for rec in alleles]
        return alleles
    
    with gzip.open(path, "r") as file:
        for trid, group in itertools.groupby(file, key=lambda line: line.decode("utf8").split()[0]):
            alleles = parse_alleles(group)
            repeat_recs = [(s, [len(a) for a in als]) for s, als in alleles]
            repeat_recs = [RepeatRec(sample, min(allele), max(allele)) for sample, allele in repeat_recs]
            
            yield trid, repeat_recs

def get_minmax_alleles(path):
    '''organize information to retrieve min and max length alleles for each loci, along with their sample IDs'''
    MinMaxRec = namedtuple("MinMaxRec", "repeat_loci min_size min_size_sample max_size max_size_sample")
    for i, (trid, repeat_recs) in enumerate(get_repeat_recs(path)):    
        loci_list = []
        min_allele = float("inf")
        max_allele = 0
        min_sample = []
        max_sample = []
        for i, rec in enumerate(repeat_recs):
            if rec.short_allele <= min_allele:
                if rec.short_allele < min_allele:
                    min_sample.clear()
                min_allele = rec.short_allele
                min_sample.append(rec.sample)
            elif rec.long_allele >= max_allele:
                if rec.long_allele > max_allele:
                    max_sample.clear()
                max_allele = rec.long_allele
                max_sample.append(rec.sample)
        yield trid, min_allele, ';'.join(min_sample), max_allele, ';'.join(max_sample)

def process_alleles():
    """
    Setting up command line argument for input file path, as well as storing all relevant min/max loci information into a pandas dataframe,
    then saving it as a compressed tsv file.
    """
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process alleles data and save to a compressed TSV file.")
    parser.add_argument(
        "input_file",
        help='''Path to the input file (e.g., /Users/Downloads/alleles.gz). The input file should contain alleles data organized in the following format:
                1) Each line of the file should represent a repeat locus
                2) The columns should include the repeat ID, sample ID, and alleles sizes separated by commas.
                Example format (per line): 'chr10_100000834_100000912_A HG00099 TTTAG,TTTAGAA'.'''
            )
    args = parser.parse_args()

    # File paths
    file_path = args.input_file

    # Incorporating date of output creation into name
    today = date.today().strftime("%Y-%m-%d")
    output_file = f"TRGT_allelesdb_{today}.tsv.gz"

    # Extract records
    records = []
    for trid, min_allele, min_sample, max_allele, max_sample in get_minmax_alleles(file_path):
        records.append([trid, min_allele, min_sample, max_allele, max_sample])

    # Create DataFrame directly from generator expression and sort based on chromosome location
    df = pd.DataFrame(records, columns=['trid', 'min_allele', 'min_sample', 'max_allele', 'max_sample']).sort_values(by='trid')

    # Save DataFrame to compressed TSV file
    with gzip.open(output_file, 'wt', encoding='utf-8') as f_out:
        df.to_csv(f_out, sep='\t', index=False, compression='gzip')

    print("TSV saved successfully.")

process_alleles()
