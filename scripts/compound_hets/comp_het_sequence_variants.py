import argparse
from functools import reduce
import pandas as pd
import numpy as np
import re

def filter_low_impact_variants(low: pd.DataFrame) -> pd.DataFrame:
    """Select low impact variants that may be damaging
    
    Args:
        low (pd.DataFrame): DataFrame with low impact variants queried from gemini database
    
    Returns:
        pd.DataFrame: DataFrame with genic low impact variants that may be damaging (SNVs with SpliceAI >=0.2, CADD >= 14, or indel with no SpliceAI or CADD score)
    """
    # filter out intergenic variants
    low_impact_var_filter = low[low['Variation'].isin(['intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'mature_miRNA_variant', 'synonymous_variant'])].copy()
    # parse SpliceAI and CADD scores
    low_impact_var_filter['SpliceAI_score_parsed'] = low_impact_var_filter["SpliceAI_score"].apply(lambda x: max(x.split("|")[2:6]) if not pd.isna(x) else np.nan)
    low_impact_var_filter['SpliceAI_score_parsed'] = low_impact_var_filter['SpliceAI_score_parsed'].astype(float)
    low_impact_var_filter['Cadd_score'] = low_impact_var_filter['Cadd_score'].astype(float)
    # identify variant type
    low_impact_var_filter["sum_ref_alt_length"] = low_impact_var_filter["Ref"].str.len() + low_impact_var_filter["Alt"].str.len()
    low_impact_var_filter["variant_type"] = low_impact_var_filter["sum_ref_alt_length"].apply(lambda x: "SNV" if x == 2 else "indel")
    # retain SNVs with SpliceAI >=0.2, CADD >= 14, or indel with no SpliceAI or CADD score:
    low_impact_var_filter_scores = low_impact_var_filter[(low_impact_var_filter["SpliceAI_score_parsed"] >= 0.2) | 
                                                        (low_impact_var_filter["Cadd_score"] >= 14) | 
                                                        ((low_impact_var_filter["variant_type"] == "indel") & (low_impact_var_filter["SpliceAI_score_parsed"].isna()) & (low_impact_var_filter["Cadd_score"].isna())) # missing spliceAI and cadd score
                                                        ]

    return low_impact_var_filter_scores

def infer_pedigree_roles(pedigree: str) -> dict:
    """Map pedigree roles to individual IDs."""
    pedigree = pd.read_csv(pedigree, sep=" ", header=None, names=["family_ID", "individual_ID", "paternal_ID", "maternal_ID", "sex", "phenotype"])
    # just take ID of first child if there are multiple children
    child = child = pedigree[(pedigree["paternal_ID"] != "0") & (pedigree["maternal_ID"] != "0")]["individual_ID"].values[0]
    father = pedigree[pedigree["individual_ID"] == pedigree[pedigree["individual_ID"] == child]["paternal_ID"].values[0]]["individual_ID"].values[0]
    mother = pedigree[pedigree["individual_ID"] == pedigree[pedigree["individual_ID"] == child]["maternal_ID"].values[0]]["individual_ID"].values[0]
    fam_dict = {"child": child, "father": father, "mother": mother}

    return fam_dict

def abstract_gt(ref, alt, gt):
    if "|" not in gt:
        return gt 
    else:
        gt_0 = gt.split("|")[0]
        gt_1 = gt.split("|")[1]
        if ref == gt_0:
            return "ref|alt"
        else:
            return "alt|ref"

def count_gene_variants_by_parental_haplotype(variant_to_gene: pd.DataFrame, variant_gt_details: pd.DataFrame, fam_dict: dict) -> pd.DataFrame: 
    """Count variants by parental haplotype for a given gene.
    
    For each heterozygous variant in the gene, determines whether it was inherited on the maternal
    or paternal haplotype based on parental genotypes. Inheritance is determined by:
    - If mother is heterozygous and father is reference, counts as maternal
    - If father is heterozygous and mother is reference, counts as paternal 
    - If one parent genotype is missing, assumes inheritance from heterozygous parent
    - If both parents are reference or missing, counts as unknown
    
    Args:
        variant_to_gene (pd.DataFrame): DataFrame mapping variant IDs to gene IDs
        variant_gt_details (pd.DataFrame): DataFrame with variant genotype details including
            zygosity for proband and parents
        fam_dict (dict): Dictionary with family roles
            
    Returns:
        pd.DataFrame: gene ids to maternal and paternal haplotype counts
    """

    gene_to_haplotype_counts = {}

    for gene in variant_to_gene["Ensembl_gene_id"].unique():
        maternal_haplotype_count = 0
        paternal_haplotype_count = 0
        unknown_haplotype_count = 0
        for variant in variant_to_gene[variant_to_gene["Ensembl_gene_id"] == gene]["Variant_id"]:
            proband_zygosity = variant_gt_details.loc[variant, fam_dict["child"]]["Zygosity"]
            mother_zygosity = variant_gt_details.loc[variant, fam_dict["mother"]]["Zygosity"]    
            father_zygosity = variant_gt_details.loc[variant, fam_dict["father"]]["Zygosity"]
            if not (proband_zygosity == "homozygous_ref" or proband_zygosity == "missing" or proband_zygosity == "homozygous_alt"):
                if not (mother_zygosity == "homozygous_alt" or father_zygosity == "homozygous_alt"):
                    if not (mother_zygosity == "heterozygous" and father_zygosity == "heterozygous"): # if both parents are heterozygous, we are probably not interested in this variant
                        if mother_zygosity == "heterozygous":
                            if father_zygosity == "homozygous_ref": # clear maternal inheritance
                                maternal_haplotype_count += 1
                            elif father_zygosity == "missing": # assume maternal inheritance
                                unknown_haplotype_count += 1
                        elif father_zygosity == "heterozygous":
                            if mother_zygosity == "homozygous_ref": # clear paternal inheritance
                                paternal_haplotype_count += 1
                            elif mother_zygosity == "missing": # assume paternal inheritance
                                unknown_haplotype_count += 1
                        elif mother_zygosity == "homozygous_ref" and father_zygosity == "homozygous_ref": # can later use phase set and phase to determine which haplotype to add to?
                            unknown_haplotype_count += 1
                        elif mother_zygosity == "missing" and father_zygosity == "missing": # can later use phase set and phase to determine which haplotype to add to?
                            unknown_haplotype_count += 1

        gene_to_haplotype_counts[gene] = [maternal_haplotype_count, paternal_haplotype_count, unknown_haplotype_count]
    gene_to_haplotype_counts = pd.DataFrame.from_dict(gene_to_haplotype_counts, orient="index")
    gene_to_haplotype_counts.columns = ["maternal_haplotype_count", "paternal_haplotype_count", "unknown_haplotype_count"]

    return gene_to_haplotype_counts

def is_compound_het(maternal_haplotype_count, paternal_haplotype_count, unknown_haplotype_count):
    """Determine if a gene is compound het based on maternal and paternal haplotype counts and unknown haplotype count."""
    if maternal_haplotype_count > 0 and paternal_haplotype_count > 0:
        return "True"
    elif unknown_haplotype_count == 1 and (maternal_haplotype_count > 0 or paternal_haplotype_count > 0):
        return "Unknown"
    elif unknown_haplotype_count > 1:
        return "Unknown"
    else:
        return "False"

def determine_compound_het_status_no_parents(variant_to_gene: pd.DataFrame, variant_gt_details: pd.DataFrame) -> pd.DataFrame: 
    """Determine compound het status for a given gene using only long-read phasing information.
    
    Args:
        variant_to_gene (pd.DataFrame): DataFrame mapping variant IDs to Ensembl gene IDs
        variant_gt_details (pd.DataFrame): DataFrame with variant genotype details
    
    Returns:
        pd.DataFrame: DataFrame with compound het status for each gene
    """
    gene_to_compound_hets = {}
    for gene in variant_to_gene["Ensembl_gene_id"].unique():
        variants = variant_gt_details[variant_gt_details["Ensembl_gene_id"] == gene]
        num_variants = len(variants) # total number of het variants in the gene 
        true_compound_hets = 0
        unknown_compound_hets = 0
        phase_sets = variants["PS"].unique()
        if num_variants > 1:
            if phase_sets.tolist() == ["."]: # variants are not phased 
                unknown_compound_hets += 1
            elif len(phase_sets) == 1: # all variants are in the same phase set
                PS = phase_sets[0]
                variants_in_PS_gts = variants[variants["PS"] == PS]["GT_abstracted"].unique()
                if "ref|alt" in variants_in_PS_gts and "alt|ref" in variants_in_PS_gts:
                    true_compound_hets += 1
            else:
                # if two phase sets, and one is unphased, can be compound het (if phased variant(s) are in the same phase set) or unknown
                # if two or more phase sets, can be compound het (if phased variant(s) are in the same phase set) or unknown
                unknown_compound_hets += 1 # if there's an unphased variant, we don't know if it's a compound het with phased variant(s)
                for PS in variants["PS"].unique():
                    if PS != ".":
                        variants_in_PS_gts = variants[variants["PS"] == PS]["GT_abstracted"].unique()
                        if "ref|alt" in variants_in_PS_gts and "alt|ref" in variants_in_PS_gts:
                            true_compound_hets += 1

        if true_compound_hets > 0:
            compound_het_status = "True"
        else:
            if unknown_compound_hets > 0:
                compound_het_status = "Unknown"
            else:
                compound_het_status = "False"
        gene_to_compound_hets[gene] = compound_het_status
    gene_to_compound_hets = pd.DataFrame.from_dict(gene_to_compound_hets, orient="index").reset_index()
    gene_to_compound_hets.columns = ["Ensembl_gene_id", "compound_het_status"]

    return gene_to_compound_hets

def melt_sample_columns(df: pd.DataFrame, id_col: str, prefix: str, value_name: str, sample_ids: list) -> pd.DataFrame:
    """Generate melted/long format dataframefor a given sample-specific field's prefix (e.g. genotype), for all samples detected"""
    sample_cols = [f"{prefix}{sample}" for sample in sample_ids]
    melted = df.melt(
        id_vars=[id_col, "Ref", "Alt"],
        value_vars=sample_cols,
        var_name="Sample",
        value_name=value_name
    )
    melted["Sample"] = melted["Sample"].str.replace(f"{prefix}", "", regex=False)

    return melted 

def get_compound_het_variants(high_impact: pd.DataFrame, variant_to_gene: pd.DataFrame, compound_het_status: pd.DataFrame, proband_id: str, sample_ids: list) -> pd.DataFrame:
    # TODO: may want to export all SNVs for cross-variant type CH analysis
    """Get all heterozygous high impactvariants for a given proband in genes annotatedwith compound het status.
    This will be used for cross-variant type CH analysis.
    
    Args:
        high_impact (pd.DataFrame): DataFrame with high impact variants
        variant_to_gene (pd.DataFrame): DataFrame mapping variant IDs to gene IDs
        compound_het_status (pd.DataFrame): DataFrame with compound het status for each gene
        proband_id (str): ID of the proband
        sample_ids (list): List of sample IDs
    
    Returns:
        pd.DataFrame: DataFrame with high impact variants that are compound het
    """
    ch_genes = compound_het_status.index.tolist()
    ch_variants = variant_to_gene[variant_to_gene["Ensembl_gene_id"].isin(ch_genes)]["Variant_id"].values.tolist()
    ch_variants_df = high_impact[high_impact["Variant_id"].isin(ch_variants)].copy()
    # add zygosity column 
    for sample in sample_ids:
        gt_type_col = f"gt_types.{sample}"
        zygosity_col = f"zygosity.{sample}"
        if gt_type_col in ch_variants_df.columns:
            ch_variants_df[zygosity_col] = ch_variants_df[gt_type_col].map(
                lambda x: "hom" if x == 3 
                        else "het" if x == 1 
                        else "-" if x == 0 
                        else "missing"
                        )
    ch_variants_df = ch_variants_df[ch_variants_df[f"zygosity.{proband_id}"] == "het"]
    # format dataframe for export
    ch_variants_df["Gnomad_af_popmax"] = ch_variants_df["Gnomad_af_popmax"].replace(-1, 0)
    ch_variants_df.drop(columns=["Variant_id", "variant_type", "sum_ref_alt_length", "GT_qual_max"], inplace=True)
    gt_type_drop = [col for col in ch_variants_df.columns if col.startswith("gt_types.")]
    gt_phases_drop = [col for col in ch_variants_df.columns if col.startswith("gt_phases.")]
    ch_variants_df.drop(columns=gt_type_drop + gt_phases_drop, inplace=True)
    ch_variants_df = ch_variants_df[["Chrom", "Pos", "Ref", "Alt"] +
     ["zygosity." + sample for sample in sample_ids] + 
     ["Variation", "Depth", "Quality", "Gene", "Clinvar", "Ensembl_gene_id", "Gnomad_af_popmax", "Cadd_score", "SpliceAI_score"] + 
     [f"gts.{sample}" for sample in sample_ids] + 
     [f"gt_alt_depths.{sample}" for sample in sample_ids] + 
     [f"gt_depths.{sample}" for sample in sample_ids] + 
     [f"gt_quals.{sample}" for sample in sample_ids] + 
     [f"PS.{sample}" for sample in sample_ids]]

    return ch_variants_df

def main():
    parser = argparse.ArgumentParser(description="Annotate compound het status for sequence variants.")
    parser.add_argument("--high_med", type=str, required=True, help="Path to HIGH-MED variants file")
    parser.add_argument("--low", type=str, required=True, help="Path to LOW variants file")
    parser.add_argument("--pedigree", type=str, required=True, help="Path to pedigree file")
    parser.add_argument("--coding_report", type=str, required=True, help="Path to coding report file")
    args = parser.parse_args()

    family = args.pedigree.split("/")[-1].split(".")[0].split("_")[0]
    high_med = pd.read_csv(args.high_med, sep='\t')
    low = pd.read_csv(args.low, sep='\t', low_memory=False)

    # process low impact variants to select high impact variants
    low_impact_var_filter_scores = filter_low_impact_variants(low)
    # create high_impact table: genic damaging noncoding variants and nonsynonymous variants
    high_impact = pd.concat([low_impact_var_filter_scores, high_med])
    # de-duplicate (may be LOW impact variants in high_med that are also in low_impact if they have ClinVar annotations)
    high_impact = high_impact.drop_duplicates(subset=["Variant_id"])
    # extract max genotype quality score for each variant 
    qual_cols = [col for col in high_impact.columns if col.startswith("gt_quals.")]
    high_impact["GT_qual_max"] = high_impact.apply(lambda x: max([x[col] for col in qual_cols]), axis=1)
    # drop unnecessary columns
    high_impact = high_impact.drop(["Old_multiallelic", "SpliceAI_score_parsed", "Nucleotide_change_ensembl", "Protein_change_ensembl"], axis=1)

    # get per-sample variant genotype status
    # extract all sample IDs based on genotype columns
    sample_ids = sorted({re.sub(r'^gts\.','', col) for col in high_impact.columns if col.startswith("gts.")})

    # split comma-separated PS string into individual PS columns for each sample
    ps_split_cols = high_impact["PS"].str.split(",", expand=True)
    for idx, sample in enumerate(sample_ids):
        high_impact[f"PS.{sample}"] = ps_split_cols[idx]
    high_impact.drop(columns=["PS"], inplace=True)

    variant_ps      = melt_sample_columns(high_impact, "Variant_id", "PS.", "PS", sample_ids)
    variant_gts     = melt_sample_columns(high_impact, "Variant_id", "gts.", "GT", sample_ids)
    variant_gt_type = melt_sample_columns(high_impact, "Variant_id", "gt_types.", "GT_type", sample_ids)
    variant_phase   = melt_sample_columns(high_impact, "Variant_id", "gt_phases.", "Phase", sample_ids)
    variant_qual    = melt_sample_columns(high_impact, "Variant_id", "gt_quals.", "GT_qual", sample_ids)

    # combine all melted dataframes on Variant_id and Sample to create one dataframe of variants in long format
    dfs = [variant_ps, variant_gts, variant_gt_type, variant_phase, variant_qual]
    # successively merges a list of DataFrames (dfs), each containing different variant/sample-level data,
    # to produce a single DataFrame (variant_gt_details) with all variant information for each (Variant_id, Sample) pair.
    variant_gt_details = reduce(lambda left, right: left.merge(right, on=["Variant_id", "Ref", "Alt", "Sample"], how="outer"), dfs)
    variant_gt_details.loc[variant_gt_details['GT'] == './.', 'GT_type'] = 2 # for some reason gemini codes some missing genotypes types as 3  
    variant_gt_details['Zygosity'] = variant_gt_details['GT_type'].map(lambda x: 'homozygous_alt' if x == 3 else 'heterozygous' if x == 1 else 'homozygous_ref' if x == 0 else 'missing') 
    variant_gt_details.set_index(["Variant_id", "Sample"], inplace=True)
    # create variant_to_gene table: variant ID and gene ID
    variant_to_gene = high_impact[["Variant_id", "Ensembl_gene_id"]].drop_duplicates()
    variant_gt_details.reset_index(inplace=True)
    variant_gt_details = variant_gt_details.merge(variant_to_gene, on="Variant_id", how="left")

    # restrict variants to those that are het in the proband
    fam_dict = infer_pedigree_roles(args.pedigree)
    proband_id = fam_dict["child"]
    variant_gt_details_proband = variant_gt_details[variant_gt_details["Sample"] == proband_id]
    variant_gt_details_proband = variant_gt_details_proband[variant_gt_details_proband["Zygosity"] == "heterozygous"]
    # abstract the genotype to be "ref|alt" or "alt|ref" to facilitate compound het status determination
    variant_gt_details_proband["GT_abstracted"] = variant_gt_details_proband.apply(lambda x: abstract_gt(x["Ref"], x["Alt"], x["GT"]), axis=1)

    # determine gene compound het status using only long-read phasing information
    compound_het_status_no_parents = determine_compound_het_status_no_parents(variant_to_gene, variant_gt_details_proband)
    print(f"CH gene status prior to parental phasing:\n {compound_het_status_no_parents.compound_het_status.value_counts()}")
    # if any gene has unknown compound het status, we attempt to determine the compound het status using parental genotypes
    # first, get genes with unknown compound het status
    unknown_gene = compound_het_status_no_parents[compound_het_status_no_parents["compound_het_status"] == "Unknown"]
    print(f"Number of genes with unknown compound het status after long-read phasing: {len(unknown_gene)}")
    unknown_variants = variant_gt_details[variant_gt_details["Ensembl_gene_id"].isin(unknown_gene["Ensembl_gene_id"])]
    unknown_variant_to_gene = variant_to_gene[variant_to_gene["Variant_id"].isin(unknown_variants["Variant_id"])]
    unknown_variants.set_index(["Variant_id", "Sample"], inplace=True)
    # use parental genotypes to determine compound het status
    gene_haplotype_counts = count_gene_variants_by_parental_haplotype(unknown_variant_to_gene, unknown_variants, fam_dict)
    gene_haplotype_counts["compound_het_status"] = gene_haplotype_counts.apply(lambda x: is_compound_het(x["maternal_haplotype_count"], x["paternal_haplotype_count"], x["unknown_haplotype_count"]), axis=1)
    print(f"CH status of unknown genes after parental phasing:\n {gene_haplotype_counts['compound_het_status'].value_counts()}")
    compound_het_status = compound_het_status_no_parents.copy().set_index("Ensembl_gene_id")
    for gene in gene_haplotype_counts.index:
        compound_het_status.loc[gene, "compound_het_status"] = gene_haplotype_counts.loc[gene, "compound_het_status"]
    compound_het_status.compound_het_status.value_counts()
    print(f"CH gene status after parental phasing:\n {compound_het_status.compound_het_status.value_counts()}")

    # export table of variants in genes with compound het status
    ch_variants_df = get_compound_het_variants(high_impact, variant_to_gene, compound_het_status, proband_id, sample_ids)
    ch_variants_df.to_csv(f"{family}_compound_het_variants.csv", index=False)

    # add compound het status to coding report
    coding_report = pd.read_csv(args.coding_report)
    compound_het_status = compound_het_status.reset_index()
    coding_report = coding_report.merge(compound_het_status, left_on="Ensembl_gene_id", right_on="Ensembl_gene_id", how="left")
    coding_report.replace({np.nan: "."}, inplace=True)
    output_filename = args.coding_report.split("/")[-1].replace(".csv", "_CH.csv")
    coding_report.to_csv(output_filename, index=False)

main()

