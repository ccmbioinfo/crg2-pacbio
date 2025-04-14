import allel
import argparse
from datetime import date
import pandas as pd
from pysam import VariantFile


def parse_sample_field(vcf_dict, vcf_df, fieldname, fieldtype):
    field_dict = {}

    for record in vcf_dict.fetch():
        for sample in record.samples:
            value = record.samples[sample][fieldname]
            value = [str(v) for v in value]
            name = record.samples[sample].name
            name = name.split(".")[0]
            if name not in field_dict.keys():
                field_dict[name] = [value]
            else:
                field_dict[name].append(value)
    for sample in record.samples:
        name = record.samples[sample].name
        name = name.split(".")[0]
        vcf_df[f"{name}_{fieldtype}"] = [
            ("|").join(field) for field in field_dict[name]
        ]

    return vcf_df


def recode_genes(disease_thresholds):
    disease_thresholds.loc[disease_thresholds["Gene"] == "DMPK/DMPKas", "Gene"] = "DMPK"
    disease_thresholds.loc[
        disease_thresholds["Gene"] == "ATXN8/ATXN8OS", "Gene"
    ] = "ATXN8"
    disease_thresholds.loc[disease_thresholds["Gene"] == "CBL2", "Gene"] = "CBL"
    disease_thresholds.loc[disease_thresholds["Gene"] == "XYLT", "Gene"] = "XYLT1"
    disease_thresholds.loc[disease_thresholds["Gene"] == "TK2/BEAN", "Gene"] = "BEAN1"
    disease_thresholds.loc[disease_thresholds["Gene"] == "MARCH6", "Gene"] = "MARCHF6"
    disease_thresholds.loc[disease_thresholds["Gene"] == "FMR1/FMR4", "Gene"] = "FMR1"
    disease_thresholds.loc[
        disease_thresholds["Gene"] == "LOC642361/NUTM2B-AS1", "Gene"
    ] = "NUTM2B-AS1"
    disease_thresholds.loc[disease_thresholds["Gene"] == "ZNF9/CNBP", "Gene"] = "CNBP" 
    disease_thresholds.loc[disease_thresholds["Gene"] == "C11orf80", "Gene"] = "C11ORF80"
    return disease_thresholds


def is_disease(motif_count, gene, threshold):
    if "_" in motif_count:
        return None
    else:
        try:
            motif_count = [int(m) for m in motif_count.split("|")]
        except:
            motif_count = "."

    is_disease = False
    if pd.isna(threshold):
        # threshold is NaN
        is_disease = None
    else:
        for count in motif_count:
            if count == ".":
                continue
            elif gene == "VWA1":
                if count == 1 or count == 3: # 2 copies is benign as per STRchive 
                    is_disease = True
            elif count >= int(threshold):
                is_disease = True
    return is_disease


def main(vcf, disease_thresholds, output_file):
    vcf_dict = allel.read_vcf(
        vcf,
        ["*"],
    )

    disease_thresholds = pd.read_csv(
        disease_thresholds,
        sep="\t",
    )

    vcf_dict["variants/ALT"] = ["|".join(alt) for alt in vcf_dict["variants/ALT"]]

    remove_keys = [
        "calldata/AP",
        "calldata/AL",
        "calldata/GT",
        "variants/FILTER_PASS",
        "variants/ID",
        "variants/altlen",
        "variants/is_snp",
        "calldata/ALLR",
        "calldata/AM",
        "calldata/MC",
        "calldata/MS",
        "calldata/SD",
        "samples",
    ]

    for key in remove_keys:
        vcf_dict.pop(key)

    vcf_df = pd.DataFrame(vcf_dict)

    # allel doesn't extract alt and ref support properly, so use pysam VariantFile instead
    vcf_dict_pysam = VariantFile(
        vcf,
    )

    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "AL", "allele_length")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "ALLR", "allele_CI")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "MC", "motif_count")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "MS", "motif_span")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "AM", "avg_methylation")
    vcf_df = parse_sample_field(vcf_dict_pysam, vcf_df, "SD", "spanning_read_support")
    vcf_df.columns = [col.replace("variants/", "") for col in vcf_df.columns]

    # recode genes names in disease threshold file for compatibility with TRGT gene names
    disease_thresholds = recode_genes(disease_thresholds)

    # extract disease threshold and disorder for each gene
    thresh_dict = {}
    disorder_dict = {}
    for index, row in disease_thresholds.iterrows():
        gene = row["Gene"]
        threshold = row["Disease threshold"]
        disorder = row["Disorder"]
        thresh_dict[gene] = threshold
        disorder_dict[gene] = disorder

    vcf_df["GENE"] = vcf_df["TRID"].str.split("_").str[1]
    vcf_df.loc[vcf_df["POS"] == 25013530, "GENE"] = "PRTS_ARX"
    vcf_df.loc[vcf_df["POS"] == 25013649, "GENE"] = "EIEE1_ARX"
    vcf_df["DISEASE_THRESHOLD"] = vcf_df["GENE"].map(thresh_dict)
    vcf_df["DISORDER"] = vcf_df["GENE"].map(disorder_dict)

    # fill in missing disease thresholds 
    vcf_df.loc[vcf_df["GENE"] == "C11ORF80", "DISEASE_THRESHOLD"] = 500.0 # OMIM
    vcf_df.loc[vcf_df["GENE"] == "TMEM185A", "DISEASE_THRESHOLD"] = 300.0 # Shaw et al 2002


    # prep report for export
    allele_length_cols = [col for col in vcf_df.columns if "allele_length" in col]
    allele_CI_cols = [col for col in vcf_df.columns if "allele_CI" in col]
    spanning_cols = [col for col in vcf_df.columns if "spanning" in col]
    motif_count_cols = [col for col in vcf_df.columns if "motif_count" in col]
    methylation_cols = [col for col in vcf_df.columns if "methylation" in col]

    # add disease outlier column
    for col in motif_count_cols:
        vcf_df[f"DISEASE_PREDICTION_{col}"] = vcf_df.apply(
            lambda row: is_disease(
                row[col],
                row["GENE"],
                row.DISEASE_THRESHOLD
            ),
            axis=1,
        )
    disease_pred_cols = [col for col in vcf_df.columns if "PREDICTION" in col]
    vcf_df["DISEASE_PREDICTION"] = vcf_df[disease_pred_cols].apply(
        lambda row: "|".join(row.values.astype(str)), axis=1
    )
    for col in disease_pred_cols:
        vcf_df.drop(col, axis=1)
    # for loci that are homozygous reference in all samples, set ALT allele to "."
    vcf_df.loc[vcf_df["ALT"] == "||", "ALT"] = "homozygous_ref"

    report_cols = (
        [
            "CHROM",
            "POS",
            "END",
            "REF",
            "ALT",
            "GENE",
            "DISORDER",
            "MOTIFS",
            "STRUC",
            "DISEASE_THRESHOLD",
            "DISEASE_PREDICTION",
        ]
        + motif_count_cols
        + allele_length_cols
        + allele_CI_cols
        + spanning_cols
        + methylation_cols
    )

    vcf_df = vcf_df[report_cols]
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    output_prefix = output_file.replace(".csv", "")
    vcf_df.to_csv(f"{output_file}", index=False)
    vcf_df.to_csv(f"{output_prefix}.{today}.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a report for known pathogenic repeat expansion loci from a merged TRGT VCF"
    )
    parser.add_argument("--vcf", type=str, help="TRGT VCF", required=True)
    parser.add_argument("--disease_thresholds", type=str, help="Repeat loci", required=True)
    parser.add_argument("--output_file", type=str, help="Output filename", required=True)

    args = parser.parse_args()
    vcf = args.vcf
    disease_thresholds = args.disease_thresholds
    output_file =args.output_file

    main(
        vcf,
        disease_thresholds,
        output_file
    )
