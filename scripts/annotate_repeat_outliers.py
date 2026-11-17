import argparse
import pandas as pd
import pyranges as pr
import numpy as np
import pysam
from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval
import sys

from annotation.annotate import (
    prepare_OMIM
)


PROMOTER_SIZE = 5000

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True, help="Path to TRGT VCF file")
    parser.add_argument("--lps", type=str, required=True, help="Path to LPS file")
    parser.add_argument("--dbscan_cutoffs_lps", type=str, required=True, help="Path to DBSCAN cutoffs for LPS")
    parser.add_argument("--dbscan_cutoffs_al", type=str, required=True, help="Path to DBSCAN cutoffs for AL")
    parser.add_argument("--hprc_stats", type=str, required=True, help="Path to HPRC stats file")
    parser.add_argument("--annotated_catalog", type=str, required=True, help="Path to TRExplorer annotated catalog file")
    parser.add_argument("--gtf", type=str, required=True, help="Path to GENCODE GTF file")
    parser.add_argument("--promoter_size", type=int, required=True, help="Promoter size")
    parser.add_argument("--SCREEN_promoters", type=str, required=True, help="Path to SCREEN promoters BED file")
    parser.add_argument("--outfile", type=str, required=True, help="Path to output file")
    parser.add_argument("--omim", type=str, required=True, help="Path to OMIM genemap2 file")
    parser.add_argument("--hpo", type=str, required=False, help="Path to HPO annotations file")
    parser.add_argument("--ped", type=str, required=True, help="Path to PED file")
    return parser.parse_args()


def is_X(trid):
    return "X" in trid or "X" in trid[0]

def process_variant(trid, variant, samples, lps_dict):
    tr_dict = {}
    tr_dict = {"TRID": trid}
    for sample in samples:
        AL = list(variant.samples[sample]['AL'])
        AM = [round(am, 2) if am is not None else None for am in variant.samples[sample]['AM']]
        if AM == [None, None]:
            AM = None
        MP = [round(mp, 2) if mp is not None else None for mp in variant.samples[sample]["AP"]]
        if MP == [None, None]:
            MP = None
        AL_hemi = None
        if is_X(trid):
            if AL is not None and len(AL) == 1: # if only one allele length reported, this is a male individual
                AL_hemi = AL[0]
        GT = list(variant.samples[sample]['GT'])
        SD = list(variant.samples[sample]['SD'])
        tr_dict[f'{sample}_LPS'] = lps_dict[trid][sample] 
        tr_dict[f'{sample}_AL'] = AL
        tr_dict[f'{sample}_AM'] = AM
        tr_dict[f'{sample}_MP'] = MP
        tr_dict[f'{sample}_AL_hemi'] = AL_hemi
        tr_dict[f'{sample}_GT'] = GT
        try:
            tr_dict[f'{sample}_DP'] = sum([int(dp) for dp in SD])
            tr_dict[f'{sample}_AD'] = SD[1] if GT != [0,0] else "."
        except:
            if is_X(trid):
                tr_dict[f'{sample}_DP'] = SD[0]
                tr_dict[f'{sample}_AD'] = SD[0] if GT != [0,0] else "."
            else:
                tr_dict[f'{sample}_DP'] = "."
                tr_dict[f'{sample}_AD'] = "."
    return tr_dict

def list_col_to_str(col):
    if col is None:
        return "."
    elif isinstance(col, int):
        return str(col)
    elif "," in col:
        return col.replace(",", "|")
    elif isinstance(col, list):
        col = [str(c) for c in col]
        return "|".join(col)
    else:
        return col

def outlier_type(AL, LPS, cutoff_AL, cutoff_LPS):
    if AL is None:
        return "False"
    
    if LPS > cutoff_LPS:
        if AL > cutoff_AL:
            return "LPS_AL_outlier"
        else:
            return "LPS_outlier"
    elif AL > cutoff_AL:
        return "AL_outlier"
    else:
        return "False"

def is_outlier(row, proband_id):

    cutoff_LPS = row["LPS_cutoff_DBSCAN"]
    cutoff_AL = row["AL_cutoff_DBSCAN"]
    cutoff_LPS_hemi = row["hemi_LPS_cutoff_DBSCAN"]
    cutoff_AL_hemi = row["hemi_AL_cutoff_DBSCAN"]
    AL = row[f'{proband_id}_AL']
    LPS = row[f'{proband_id}_LPS']
    
    # handle hemizygous genotypes
    if len(AL.split("|")) == 1 and row["TRID"].startswith("X"): 
        if AL != "None" and AL != ".":
            AL = int(AL)
        else:
            AL = 0
        if LPS != "None" and LPS != ".":
            LPS = int(LPS)
        else:
            LPS = 0
        return outlier_type(AL, LPS, cutoff_AL_hemi, cutoff_LPS_hemi)
    else:
        try:
            AL = sum([int(al) for al in AL.split("|") if al is not None])
        except:
            AL = None
        lps = row[f'{proband_id}_LPS']
        if lps != ".":
            LPS = sum([int(lps) for lps in LPS.split("|")])
        else:
            LPS = 0

        return outlier_type(AL, LPS, cutoff_AL, cutoff_LPS)

def calc_z_score(x, mean, std):
    # get max of the two alleles
    if x != "." and not pd.isna(x) and x != "None":
        x = max([int(lps) for lps in x.split("|")])
        mean = float(mean)
        std = float(std)
        try:
            return round((x - mean) / std, 3)
        except ZeroDivisionError:
            return round((x - mean) / (std + 0.1), 3)
        except TypeError:
            print(f"TypeError: {x} - {mean} / {std}")
            return np.nan
    else:
        return np.nan

def get_distance_from_cutoff(row, proband_id, stat):
    trid = row["TRID"]
    if stat == "LPS":
        if row[f'{proband_id}_LPS'] == ".":
            return np.nan
        LPS = sum([int(lps) for lps in row[f'{proband_id}_LPS'].split("|")])
        if "X" not in trid:
            return round(LPS - float(row["LPS_cutoff_DBSCAN"]), 3)
        else:
            return round(LPS - float(row["hemi_LPS_cutoff_DBSCAN"]), 3)
    elif stat == "AL":
        if row[f'{proband_id}_AL'] == 'None':
            return np.nan
        AL = sum([int(al) for al in row[f'{proband_id}_AL'].split("|")])
        if "X" not in trid:
            return round(AL - float(row["AL_cutoff_DBSCAN"]), 3)
        else:
            return round(AL - float(row["hemi_AL_cutoff_DBSCAN"]), 3)

def is_greater_than_max_HPRC_AoU1027(row, id, outlier_type):
    TRID = row["TRID"]
    if outlier_type == "False":
        return "False"
    if outlier_type == "LPS_outlier" or outlier_type == "LPS_AL_outlier":
        LPS = row[f'{id}_LPS'].split("|")
        if "X" in TRID and len(LPS) == 1:
            LPS = int(LPS[0])
            lps_col = "lps_max_hemi"
        else:
            LPS = max([int(lps) for lps in LPS])
            lps_col = "lps_max"
        try:
            HPRC_max = int(row[lps_col])
        except:
            HPRC_max = 0
        if LPS > HPRC_max:
            AoU1027_lps_max = row["AoU1027_MaxAllele"] # we only have max LPS for AoU1027
            if not pd.isna(AoU1027_lps_max):
                if LPS > AoU1027_lps_max:
                    return "True"
                else:
                    return "False"
            else:
                return "True"
        else:
            return "False"
    else:
        AL = row[f'{id}_AL'].split("|")
        if "X" in TRID and len(AL) == 1:
            AL = int(AL[0])
            al_col = "allele_len_max_hemi"
        else:
            AL = max([int(al) for al in AL])
            al_col = "allele_len_max"
        try:
            HPRC_max = int(row[al_col])
        except:
            HPRC_max = 0
        if AL > HPRC_max:
            return "True"
        else:
            return "False"

def compute_genomic_region_of_row(row, genes_gtf_path=None, verbose=False, show_progress_bar=False, n=None, promoter_size=PROMOTER_SIZE):
    TRID = row["TRID"]
    chrom = "chr" + str(TRID.split("-")[0])
    start_0based = TRID.split("-")[1]
    end = TRID.split("-")[2]
    
    return compute_genomic_region_of_interval(chrom, int(start_0based) + 1, int(end), genes_gtf_path=genes_gtf_path, verbose=verbose, show_progress_bar=show_progress_bar, n=n, promoter_size=promoter_size)

def main():
    args = parse_args()
    # Create a dictionary of TRID to LPS values for each sample
    print("Creating dictionary of TRID to LPS values for each sample")
    lps = pd.read_csv(args.lps, sep="\t", engine="pyarrow")
    lps = lps[~lps["trid"].str.contains("VC:")]
    samples = [col for col in lps.columns if col != "trid" and col != "motif"]
    lps_dict = lps.set_index("trid")[samples].to_dict("index")

    # Parse VCF file to extract allele lengths and other information
    print("Parsing VCF file to extract allele lengths and other information")
    with pysam.VariantFile(args.vcf, 'r') as vcf_in:
        samples = list(vcf_in.header.samples)
        variant_df = pd.DataFrame()
        for variant in vcf_in:
            trid = variant.info['TRID']
            if not "VC" in trid:
                allele_dict = process_variant(trid, variant, samples, lps_dict)
                variant_df = pd.concat([variant_df, pd.DataFrame([allele_dict])], ignore_index=True)

    # Convert list columns into "|" separated strings
    cols = [col for col in variant_df.columns if col != "TRID"]
    for col in cols:
        variant_df[col] = variant_df[col].apply(lambda x: list_col_to_str(x))

    # Filter out variation clusters
    variant_df = variant_df[~variant_df["TRID"].str.contains("VC:")]

    # Import DBSCAN cutoffs
    LPS_cutoffs = pd.read_csv(args.dbscan_cutoffs_lps, sep="\t", engine="pyarrow")
    LPS_cutoffs.columns = ["TRID", "LPS_cutoff_DBSCAN"]
    LPS_cutoffs_hemi = LPS_cutoffs[LPS_cutoffs["TRID"].str.startswith("hemi")].rename(columns={"LPS_cutoff_DBSCAN": "hemi_LPS_cutoff_DBSCAN"})
    LPS_cutoffs_hemi["TRID"] = LPS_cutoffs_hemi["TRID"].str.replace("hemi", "")
    AL_cutoffs = pd.read_csv(args.dbscan_cutoffs_al, sep="\t", engine="pyarrow")
    AL_cutoffs.columns = ["TRID", "AL_cutoff_DBSCAN"]
    AL_cutoffs_hemi = AL_cutoffs[AL_cutoffs["TRID"].str.startswith("hemi")].rename(columns={"AL_cutoff_DBSCAN": "hemi_AL_cutoff_DBSCAN"})
    AL_cutoffs_hemi["TRID"] = AL_cutoffs_hemi["TRID"].str.replace("hemi", "")

    # Merge TRs with DBSCAN cutoffs
    variant_df_cutoffs = variant_df.merge(LPS_cutoffs, on="TRID", how="left").merge(AL_cutoffs, on="TRID", how="left")
    # Add hemizygous cutoffs
    variant_df_cutoffs = variant_df_cutoffs.merge(LPS_cutoffs_hemi, on="TRID", how="left").merge(AL_cutoffs_hemi, on="TRID", how="left")

    # Filter out non-variant calls in proband
    proband_id = "03"
    variant_df_cutoffs = variant_df_cutoffs[(variant_df_cutoffs[f"{proband_id}_GT"] != "0|0") & (variant_df_cutoffs[f"{proband_id}_GT"] != "None")].copy()
    # Filter out loci not present in cutoffs file, i.e. not in HPRC v2.0 VCF
    variant_df_cutoffs = variant_df_cutoffs[variant_df_cutoffs["AL_cutoff_DBSCAN"].notna()]

    # Find outliers
    print("Finding outliers")
    for id in samples:
        variant_df_cutoffs[f"{id}_OUTLIER"] = variant_df_cutoffs.apply(lambda row: is_outlier(row, id), axis=1)

    # Add HPRC stats
    print("Adding HPRC stats")
    stats = pd.read_csv(args.hprc_stats, sep="\t", compression="gzip", engine="pyarrow")
    variant_df_cutoffs = variant_df_cutoffs.merge(stats, left_on="TRID", right_on="trid", how="left")

    variant_df_cutoffs = variant_df_cutoffs[variant_df_cutoffs["AN"].notna()] # TODO: exclude this filter for cohort analyses?


    for id in samples:
        variant_df_cutoffs[f'{id}_distance_from_LPS_cutoff_DBSCAN'] = variant_df_cutoffs.apply(lambda row: get_distance_from_cutoff(row, id, "LPS"), axis=1)
        variant_df_cutoffs[f'{id}_distance_from_AL_cutoff_DBSCAN'] = variant_df_cutoffs.apply(lambda row: get_distance_from_cutoff(row, id, "AL"), axis=1)
        variant_df_cutoffs[f'{id}_z_score_LPS'] = variant_df_cutoffs.apply(lambda row: calc_z_score(row[f'{id}_LPS'], row["lps_mean"], row["lps_std"]), axis=1)
        variant_df_cutoffs[f'{id}_z_score_AL'] = variant_df_cutoffs.apply(lambda row: calc_z_score(row[f'{id}_AL'], row["allele_len_mean"], row["allele_len_std"]), axis=1)
        variant_df_cutoffs[f'{id}_z_score_max'] = variant_df_cutoffs[[f'{id}_z_score_LPS', f'{id}_z_score_AL']].max(axis=1)

    # Add motif column
    variant_df_cutoffs["MOTIF"] = variant_df_cutoffs["TRID"].str.split("-").str[3]


    # Add TRExplorer annotations and stats
    print("Adding TRExplorer annotations and stats")
    anno_cols = ['LocusId',
        'KnownDiseaseAssociatedMotif', 'KnownDiseaseAssociatedLocus',
        'NumRepeatsInReference',
        'ReferenceRepeatPurity',
        'AoU1027_MaxAllele', 'AoU1027_Stdev',
        'AoU1027_99thPercentile', 'AoU1027_NumCalledAlleles',
        'AoU1027_OE_Length','AoU1027_OE_LengthPercentile',
        'TenK10K_MaxAllele', 'TenK10K_Stdev', 'TenK10K_99thPercentile', 'TenK10K_NumCalledAlleles']

    annos = pd.read_csv(args.annotated_catalog, sep="\t", compression="gzip", usecols=anno_cols, engine="pyarrow")
    variant_df_cutoffs_anno = variant_df_cutoffs.merge(annos, left_on="TRID", right_on="LocusId", how="left")

    # Add GENCODE gene annotations
    print("Adding GENCODE gene annotations")
    variant_df_cutoffs_anno[[f"GencodeGeneRegion", f"GencodeGeneName", f"GencodeGeneId", f"GencodeTranscriptId"]] = variant_df_cutoffs_anno.apply(
        lambda input_row, gtf=args.gtf: compute_genomic_region_of_row(input_row, genes_gtf_path=gtf, verbose=False, show_progress_bar=False, n=None, promoter_size=PROMOTER_SIZE), axis=1, result_type='expand')

    # Add OMIM phenotype
    print("Adding OMIM phenotype")
    omim = prepare_OMIM(f"{args.omim}/genemap2.txt")
    omim = omim[omim["gene_id"].notna()]
    variant_df_cutoffs_anno = variant_df_cutoffs_anno.merge(omim, how="left", left_on="GencodeGeneId", right_on="gene_id")

    # Add HPO annotations
    print("Adding HPO annotations")
    if args.hpo is not None:
        hpo = pd.read_csv(args.hpo, sep="\t", usecols=["Gene ID", "Features"])
        hpo = hpo[hpo["Gene ID"].notna()]
        variant_df_cutoffs_anno = variant_df_cutoffs_anno.merge(hpo, how="left", left_on="GencodeGeneId", right_on="Gene ID")
    else:
        variant_df_cutoffs_anno["Features"] = "."

    # Add SCREEN promoters
    print("Adding SCREEN promoters")
    promoters = pr.read_bed(args.SCREEN_promoters)
    # Convert STRs to BED format and then PyRanges object
    variant_bed = variant_df_cutoffs_anno[["TRID"]].copy()
    variant_bed["Start"] = variant_bed["TRID"].str.split("-").str[1]
    variant_bed["End"] = variant_bed["TRID"].str.split("-").str[2]
    variant_bed["Chromosome"] = variant_bed["TRID"].str.split("-").str[0]
    variant_pr = pr.PyRanges(variant_bed)
    # Intersect STRs with promoters
    variant_promoters = variant_pr.join(promoters, how="left", suffix="_promoter").df
    variant_promoters = variant_promoters.rename(columns={"Name": "SCREEN_ID", "Score": "SCREEN_promoter_ID", "Strand": "Type"})
    variant_promoters = variant_promoters[["TRID", "SCREEN_promoter_ID"]].replace("-1", ".")
    # Merge back to variant dataframe
    variant_df_cutoffs_anno = variant_df_cutoffs_anno.merge(variant_promoters, on="TRID", how="left")

    # Is AL or LPS greater than max HPRC or AoU1027?
    print("Is AL or LPS greater than max HPRC or AoU1027?")
    for id in samples:
        variant_df_cutoffs_anno[f'{id}_greater_than_max_HPRC_AoU1027'] = variant_df_cutoffs_anno.apply(lambda row: is_greater_than_max_HPRC_AoU1027(row, id, row[f"{id}_OUTLIER"]), axis=1)

    # Format for export
    print("Formatting for export")
    variant_df_cutoffs_anno.rename(columns={"AN": "HPRC_AN", 'allele_len_mean': "HPRC_allele_len_mean",
        'allele_len_std': "HPRC_allele_len_std", 'allele_len_99': "HPRC_allele_len_99", 'allele_len_max': "HPRC_allele_len_max",
        'allele_len_hemi_99': "HPRC_allele_len_hemi_99", 'allele_len_max_hemi': "HPRC_allele_len_max_hemi", 'AM_mean': "HPRC_AM_mean", 'AM_std': "HPRC_AM_std",
        'MP_mean': "HPRC_MP_mean", 'MP_std': "HPRC_MP_std", 'lps_mean': "HPRC_lps_mean", 'lps_std': "HPRC_lps_std", 'lps_99': "HPRC_lps_99", 'lps_max': "HPRC_lps_max",
        'lps_hemi_99': "HPRC_lps_hemi_99", 'lps_max_hemi': "HPRC_lps_max_hemi",
        'AoU1027_MaxAllele': "AoU1027_lps_max", 'AoU1027_Stdev': "AoU1027_lps_std", 'AoU1027_99thPercentile': "AoU1027_lps_99",
        'AoU1027_NumCalledAlleles': "AoU1027_AN", 'TenK10K_MaxAllele': "TenK10K_lps_max",
        'TenK10K_NumCalledAlleles': "TenK10K_AN", 'TenK10K_Stdev': "TenK10K_lps_std", 'TenK10K_99thPercentile': "TenK10K_lps_99" }, inplace=True)

    samples_cols = []
    for id in samples:
        samples_cols.append(f"{id}_OUTLIER")
        samples_cols.append(f"{id}_z_score_max")
        samples_cols.append(f"{id}_distance_from_LPS_cutoff_DBSCAN")
        samples_cols.append(f"{id}_distance_from_AL_cutoff_DBSCAN")
        samples_cols.append(f"{id}_GT")
        samples_cols.append(f"{id}_LPS")
        samples_cols.append(f"{id}_AL")
        samples_cols.append(f"{id}_AM")
        samples_cols.append(f"{id}_MP")
        samples_cols.append(f"{id}_AD")
        samples_cols.append(f"{id}_DP")
        samples_cols.append(f"{id}_greater_than_max_HPRC_AoU1027")

    cutoff_cols = ["LPS_cutoff_DBSCAN", "AL_cutoff_DBSCAN", "hemi_LPS_cutoff_DBSCAN", "hemi_AL_cutoff_DBSCAN"]

    variant_df_cutoffs_anno = variant_df_cutoffs_anno[['TRID','MOTIF'] + ['GencodeGeneName', 'GencodeGeneRegion', 'GencodeGeneId',  'omim_phenotype', 'omim_inheritance', 'Features','SCREEN_promoter_ID', 'KnownDiseaseAssociatedMotif', 'KnownDiseaseAssociatedLocus'] + cutoff_cols + samples_cols +  ['HPRC_AN', 'HPRC_allele_len_mean',
        'HPRC_allele_len_std', 'HPRC_allele_len_99', 'HPRC_allele_len_max',
        'HPRC_allele_len_hemi_99', 'HPRC_allele_len_max_hemi', 'HPRC_AM_mean',
        'HPRC_AM_std', 'HPRC_MP_mean', 'HPRC_MP_std', 'HPRC_lps_mean',
        'HPRC_lps_std', 'HPRC_lps_99', 'HPRC_lps_max', 'HPRC_lps_hemi_99',
        'HPRC_lps_max_hemi',
        'AoU1027_lps_max', 'AoU1027_lps_std', 'AoU1027_lps_99',
        'AoU1027_AN', 'AoU1027_OE_Length','AoU1027_OE_LengthPercentile', 'TenK10K_lps_max',
        'TenK10K_AN', 'TenK10K_lps_std', 'TenK10K_lps_99' ]]

    variant_df_cutoffs_anno.rename(columns={"omim_phenotype": "OMIM_phenotype", "omim_inheritance": "OMIM_inheritance", "Features": "HPO",
    "GencodeGeneName": "GENCODE_gene_name", "GencodeGeneRegion": "GENCODE_gene_region", "GencodeGeneId": "GENCODE_gene_id",
    "KnownDiseaseAssociatedMotif": "Known_disease_asociated_motif", "KnownDiseaseAssociatedLocus": "Known_disease_associated_locus"}, inplace=True)

    # Filter for outliers in proband
    variant_df_cutoffs_anno_outliers = variant_df_cutoffs_anno[variant_df_cutoffs_anno[f"{proband_id}_OUTLIER"] != "False"].copy()

    # Filter out variants with allele length less than max HPRC
    variant_df_cutoffs_anno_outliers = variant_df_cutoffs_anno_outliers[variant_df_cutoffs_anno_outliers[f"{proband_id}_greater_than_max_HPRC_AoU1027"] == "True"]

    variant_df_cutoffs_anno_outliers.replace(np.nan, ".", inplace=True)
    variant_df_cutoffs_anno_outliers = variant_df_cutoffs_anno_outliers.sort_values(by=[f"{proband_id}_z_score_max"], ascending=False)

    # Filter out intergenic variants
    variant_df_cutoffs_anno_outliers = variant_df_cutoffs_anno_outliers[variant_df_cutoffs_anno_outliers["GENCODE_gene_region"] != "intergenic"]
    

    # Write to file 
    print(f"Writing to file {args.outfile}")
    variant_df_cutoffs_anno_outliers.to_csv(args.outfile, index=False)

if __name__ == "__main__":
    main()