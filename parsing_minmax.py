##Developing a Database that tells us the min&max repeat sizes in samples for each locus
#Columns: repeat	min_size	max_size	min_size_sample	max_size_sample (5 columns)

##Steps
#Iterate over each line and group everything by chromosome loci
#for each allele (2 of them per sample), print out the length
#For each loci, extract the min and max allele lengths, along with their sampleIDs (5 values, loci, min, max, sampleIDmin, sampleIDmax)
#Create a dataframe of this information (pandas) and sort them based on chromosome location (ascending order)


import gzip
import itertools
import numpy as np
import pandas as pd
from collections import namedtuple

RepeatRec = namedtuple("RepeatRec", "sample short_allele long_allele")

#Parse through file and extract the repeat information (sample, alleles' size) along with their chrom location
def get_repeat_recs(path):
    def parse_alleles(group):
        alleles = list(line.decode("utf8").split() for line in group)
        alleles = [(rec[1], rec[2].split(",")) for rec in alleles]
        return alleles
    
    with gzip.open(path, "r") as file:
        for trid, group in itertools.groupby(file, key=lambda line: line.decode("utf8").split()[0]):
            alleles = parse_alleles(group)
            repeat_recs = [(s, [len(a) for a in als]) for s, als in alleles]
            repeat_recs = [RepeatRec(s, min(als), max(als)) for s, als in repeat_recs]
            
            yield trid, repeat_recs

MinMaxRec = namedtuple("MinMaxRec", "repeat_loci min_size min_size_sample max_size max_size_sample")

#organize information to retrieve min and max length alleles for each loci, along with their sample IDs
def get_minmax_alleles(path):
    for i, (trid, repeat_recs) in enumerate(get_repeat_recs(path)):    
        min_allele = float("inf")
        max_allele = 0
        min_sample = ""
        max_sample = ""
        for i, rec in enumerate(repeat_recs):
            if rec.short_allele < min_allele or rec.short_allele == min_allele:
                min_allele = rec.short_allele
                min_sample = rec.sample
            elif rec.long_allele > max_allele:
                max_allele = rec.long_allele
                max_sample = rec.sample
        yield trid, min_allele, min_sample, max_allele, max_sample

allele_minmax_db = []
for trid, min_allele, min_sample, max_allele, max_sample in get_minmax_alleles("../Downloads/alleles200.gz"):
    allele_minmax_db.append(MinMaxRec(trid, min_allele, min_sample, max_allele, max_sample))


#Create a dataframe of this information (pandas) and sort them based on chromosome location 
df = pd.DataFrame(allele_minmax_db)

# Save DataFrame to compressed TSV file
df.to_csv('TRGT_alleles_db.tsv.gz', sep='\t', index=False, compression='gzip')