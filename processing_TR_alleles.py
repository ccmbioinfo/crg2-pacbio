import pandas as pd
from datetime import date
import argparse 

def process_allele(input_filepath):
    '''
    Importing a file containing tandem repeat (TR) alleles and parsing through them to retrieve the smallest and largest 
    repeat sizes within each genomic loci, along with the sampleIDs that contain these repeats.
    '''
    #Processing input data into desired output
    alleles_df = pd.read_csv(input_filepath, sep= ' ', header=None, names=["TR_ID", "SAMPLE_ID", "TR_ALLELES"])

    #Creating columns for allele_length
    alleles_df[["TR_ALLELE1", "TR_ALLELE2"]] = alleles_df["TR_ALLELES"].str.split(',', expand=True)
    alleles_df["TR_ALLELE1_LEN"] = alleles_df["TR_ALLELE1"].apply(len)
    alleles_df["TR_ALLELE2_LEN"] = alleles_df["TR_ALLELE2"].apply(len)

    #Extracting min/max alleles and their associated sampleIDs and storing into dataframe
    trid_list = []
    for trid, sample_repeats in alleles_df.groupby("TR_ID"):
        min_allele = min(min(sample_repeats["TR_ALLELE1_LEN"]), min(sample_repeats["TR_ALLELE2_LEN"]))
        max_allele = max(max(sample_repeats["TR_ALLELE1_LEN"]), max(sample_repeats["TR_ALLELE2_LEN"]))
        min_samples = ';'.join(sample_repeats.loc[(sample_repeats["TR_ALLELE1_LEN"] == min_allele) | (sample_repeats["TR_ALLELE2_LEN"] == min_allele),'SAMPLE_ID'])
        max_samples = ';'.join(sample_repeats.loc[(sample_repeats["TR_ALLELE1_LEN"] == max_allele) | (sample_repeats["TR_ALLELE2_LEN"] == max_allele),'SAMPLE_ID'])
        trid_list.append({
            "TR_ID": trid,
            "MIN_ALLELE_SIZE" : min_allele,
            "MIN_SAMPLE_ID" : min_samples,
            "MAX_ALLELE_SIZE" : max_allele,
            "MAX_SAMPLE_ID" : max_samples
        })
    
    final_df = pd.DataFrame(trid_list)
    return final_df
    
def input_output_handling():
    # Setting up command-line argument for input_path
    parser = argparse.ArgumentParser(description="Process repeat allele data and save to a g-zipped TSV file.")
    parser.add_argument(
        "input_filepath",
        help='''Path to the input file (e.g., /Users/Downloads/alleles.gz). The input file should contain alleles data organized in the following format:
                1) Each line of the file should represent a repeat locus
                2) The columns should include the repeat ID, sample ID, and alleles sizes separated by commas.
                Example format (per line): 'chr10_100000834_100000912_A HG00099 TTTAG,TTTAGAA'.'''
            )
    args = parser.parse_args()
    input_filepath = args.input_filepath

    #Establishing output file
    today = date.today().strftime("%Y-%m-%d")
    output_file = f'TR_allele_size_db_{today}.gz'

    final_df = process_allele(input_filepath)

    try:
        final_df.to_csv(output_file, sep='\t', index=False, compression='gzip')
        print("TSV saved successfully.")
    except Exception as e:
        print(f"Error occurred while saving the TSV file: {str(e)}")

input_output_handling()