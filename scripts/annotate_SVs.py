import re
import pandas as pd
import numpy as np
from pysam import VariantFile
import argparse
from annotation.annotate import prepare_OMIM, annotate_OMIM
from collections import defaultdict
from pybedtools import BedTool
from sigfig import round
from datetime import date
import io
import sys

def rename_SV_cols(annotsv_df):
    annotsv_df.rename(
        columns={
            "SV_chrom": "CHROM",
            "SV_start": "POS",
            "SV_end": "END",
            "SV_length": "SVLEN",
            "SV_type": "SVTYPE",
        },
        inplace=True,
    )
    return annotsv_df


def annotsv_df_to_bed(annotsv_df):
    annotsv_df = annotsv_df[["CHROM", "POS", "END", "SVTYPE", "ID"]]
    bed = BedTool.from_dataframe(annotsv_df)
    return bed


def filter_length(row):
    length = row.SVLEN
    if pd.isna(length):
        return True
    elif abs(length) >= 50:
        return True
    else:
        return False


def apply_filter_length(annotsv_df):
    annotsv_df_filtered = annotsv_df.apply(lambda row: filter_length(row), axis=1)
    return annotsv_df_filtered


def filter_benign(row):
    svtype = row.SVTYPE
    loss_AF_max = row.B_loss_AFmax
    gain_AF_max = row.B_gain_AFmax
    ins_AF_max = row.B_ins_AFmax
    inv_AF_max = row.B_inv_AFmax
    if svtype == "DEL":
        if loss_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "DUP":
        if gain_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "INS":
        if ins_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "INV":
        if inv_AF_max == "nan":
            return True
        else:
            return False
    # no benign annotations for BNDs
    else:
        return True


def apply_filter_benign(annotsv_df):
    annotsv_df_filtered = annotsv_df.apply(lambda row: filter_benign(row), axis=1)
    return annotsv_df_filtered


def merge_full_split_annos(annotsv_df):
    annotsv_split = annotsv_df[annotsv_df["Annotation_mode"] == "split"]
    annotsv_full = annotsv_df[annotsv_df["Annotation_mode"] == "full"]
    DDD_cols = ["DDD_status", "DDD_mode", "DDD_consequence", "DDD_disease", "DDD_pmid"]
    GenCC_cols = ["GenCC_disease", "GenCC_moi", "GenCC_classification", "GenCC_pmid"]
    transcript_cols = [
        "Tx",
        "Tx_start",
        "Tx_end",
        "Overlapped_tx_length",
        "Overlapped_CDS_length",
        "Overlapped_CDS_percent",
        "Frameshift",
        "Exon_count",
        "Location",
        "Location2",
        "Dist_nearest_SS",
        "Nearest_SS_type",
        "Intersect_start",
        "Intersect_end",
    ]
    cols_of_interest = DDD_cols + GenCC_cols + transcript_cols
    annotsv_split_agg = (
        annotsv_split.groupby("AnnotSV_ID")
        .agg(
            {
                "GenCC_disease": ";".join,
                "GenCC_moi": ";".join,
                "GenCC_classification": ";".join,
                "GenCC_pmid": ";".join,
                "Tx": ";".join,
                "Tx_start": ";".join,
                "Tx_end": ";".join,
                "Overlapped_tx_length": ";".join,
                "Overlapped_CDS_length": ";".join,
                "Overlapped_CDS_percent": ";".join,
                "Frameshift": ";".join,
                "Exon_count": ";".join,
                "Location": ";".join,
                "Location2": ";".join,
                "Nearest_SS_type": ";".join,
                "Dist_nearest_SS": ";".join,
                "Intersect_start": ";".join,
                "Intersect_end": ";".join,
                "DDD_status": ";".join,
                "DDD_mode": ";".join,
                "DDD_consequence": ";".join,
                "DDD_disease": ";".join,
                "DDD_pmid": ";".join,
            }
        )
        .reset_index()
    )

    annotsv_full = annotsv_full[
        [col for col in annotsv_full.columns if col not in cols_of_interest]
    ]
    annotsv_merge = annotsv_split_agg.merge(
        annotsv_full, how="outer", on="AnnotSV_ID"
    ).reset_index()
    return annotsv_merge


def add_hpo(hpo, gene):
    try:
        genes = gene.split(";")
    except AttributeError:
        return "NA"
    terms = []
    for gene in genes:
        # split by - for intergenic variants, which are annotated as <upstream_gene>-<downstream_gene>
        gene = gene.split("-")
        for g in gene:
            try:
                term = str(hpo[hpo["Gene ID"] == g]["Features"].values[0])
                term = term.replace("; ", ";").split(";")
                term = list(set(term))
                for t in term:
                    if t not in terms:
                        terms.append(t)
            except IndexError:
                pass
    if len(terms) == 0:
        return "nan"
    else:
        terms = ",".join(terms)
        return terms


def add_omim(omim_df, gene):
    genes = gene.split(";")
    phenos = []
    inheritance = []
    for gene in genes:
        # split by - for intergenic variants, which are annotated as <upstream_gene>-<downstream_gene>
        gene = gene.split("-")
        for g in gene:
            try:
                phenos.append(
                    str(omim_df[omim_df["gene_id"] == g]["omim_phenotype"].values[0])
                )
                inheritance.append(
                    str(
                        omim_df[omim_df["gene_id"] == g]["omim_inheritance"].values[0]
                    )
                )
            except IndexError:
                pass
    if len(phenos) != 0:
        phenos = ",".join(phenos)
    else:
        phenos = "."
    if len(inheritance) != 0:
        inheritance = ",".join(inheritance)
    else:
        inheritance = "."
    return [phenos, inheritance]


def get_genotype(sample_GT_AD_DP):
    genotype = sample_GT_AD_DP.split(":")[0]
    if genotype == "0/1" or genotype == "0|1" or genotype == "1|0":
        zyg = "het"
    elif genotype == "1/1":
        zyg =  "hom"
        genotype = "1//1"
    elif genotype == "./.":
        zyg = genotype
    elif genotype == "0/0":
        zyg =  "-"
    else:
        zyg = genotype
    return [zyg, genotype]

def get_phase_set(sample_GT_AD_DP):
    try:
        ps = sample_GT_AD_DP.split(":")[3]
    except IndexError:
        ps = "."
    return ps

def get_alt_depth(sample_GT_AD_DP):
    try:
        alt_depth = sample_GT_AD_DP.split(":")[1].split(",")[1]
    except:
        alt_depth = 0
    return alt_depth


def get_depth(sample_GT_AD_DP):
    depth = sample_GT_AD_DP.split(":")[2]
    return depth


def get_exon_counts(annotsv_df, exon_bed):
    # original functions from SVRecords/SVAnnotator.py
    # this function includes the PB SV ID, because there can be multiple SVs with the same CHROM POS END SVTYPE that have different ALTs
    # these are uniquely identified by the PB SV ID
    exon_counts = defaultdict(
        int
    )  # incrementing dict and setting value in df is faster than incrementing values in the df
    exon_ref = BedTool(exon_bed)
    sample_bedtool = BedTool(
        list(annotsv_df.reset_index()[["CHROM", "POS", "END", "SVTYPE", "ID"]].values)
    )

    for interval in sample_bedtool.intersect(exon_ref, wa=True):
        exon_counts[
            (
                str(interval.chrom),
                str(interval.start),
                str(interval.stop),
                str(interval[3]),
                str(interval[4]),
            )
        ] += 1

    count_df = pd.Series(exon_counts).to_frame()
    count_df = count_df.reset_index()
    count_df.columns = ["CHROM", "POS", "END", "SVTYPE", "ID", "EXONS_SPANNED"]
    count_df = count_df.astype(str)
    annotsv_df = pd.merge(
        annotsv_df, count_df, how="left", on=["CHROM", "POS", "END", "SVTYPE", "ID"]
    ).fillna(value={"EXONS_SPANNED": 0})
    return annotsv_df


def annotate_pop_svs(annotsv_df, pop_svs, cols):
    # intersect annotsv and population SV bed file
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    pop_svs = pd.read_csv(pop_svs, sep="\t")
    # separate out INS and BND calls for window-based matching, DEL, DUP, and INV calls for SVLEN-based matching
    pop_svs_ins_bnd = pop_svs[(pop_svs['SVTYPE'] == 'INS') | (pop_svs['SVTYPE'] == 'BND')]
    pop_svs_del_dup_inv = pop_svs[(pop_svs['SVTYPE'] == 'DEL') | (pop_svs['SVTYPE'] == 'DUP') | (pop_svs['SVTYPE'] == 'INV')]

    # look for population INS or BND calls within 10 bp of sample INS or BND calls
    pop_svs_ins_bnd_bed = BedTool.from_dataframe(pop_svs_ins_bnd)
    window_ins_bnd = annotsv_bed.window(pop_svs_ins_bnd_bed, w=10).to_dataframe()
    window_ins_bnd.columns = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "ID",
        "CHROM_pop",
        "POS_pop",
        "END_pop",
        "SVTYPE_pop",
    ] + cols

   # look for population DEL, DUP, and INV calls with 50% reciprocal overlap with sample DEL, DUP, and INV calls
    pop_svs_del_dup_inv_bed = BedTool.from_dataframe(pop_svs_del_dup_inv)
    intersect = annotsv_bed.intersect(
        pop_svs_del_dup_inv_bed, wa=True, wb=True, F=0.5, f=0.5
    ).to_dataframe()
    intersect.columns = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "ID",
        "CHROM_pop",
        "POS_pop",
        "END_pop",
        "SVTYPE_pop",
    ] + cols


    # now concatenate the window and SVLEN-based matches
    intersect = pd.concat([window_ins_bnd, intersect])
    # popSV and sample SV must be same type
    intersect = intersect[intersect["SVTYPE"] == intersect["SVTYPE_pop"]]

    # make a column with SV details, e.g DEL:1:25266309-25324509
    pop_name = cols[0].split("_")[0]
    intersect[f"{pop_name}_SV"] = intersect[
        ["CHROM_pop", "POS_pop", "END_pop", "SVTYPE_pop"]
    ].apply(lambda x: f"{x[3]}:{x[0]}:{x[1]}-{x[2]}", axis=1)
    cols.append(f"{pop_name}_SV")
    intersect = intersect[["CHROM", "POS", "END", "SVTYPE", "ID"] + cols]
    # round AFs
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[AF_col] = [
            round(float(af), sigfigs=3) for af in intersect[AF_col].values
        ]
    except IndexError:
        pass
    # group by SV, joining annotation columns
    intersect = intersect.astype(str)
    intersect = (
        intersect.groupby(["CHROM", "POS", "END", "SVTYPE", "ID"])[cols]
        .agg({col: "; ".join for col in cols})
        .reset_index()
    )
    # get max allele frequency
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[f"{pop_name}_maxAF"] = intersect[AF_col].apply(
            lambda x: max([float(af) for af in x.split("; ")])
        )
        cols.append(f"{pop_name}_maxAF")
    # get max allele counts for C4R
    except IndexError:
        try:
            count_cols = ["C4R_AC", "seen_in_C4R_count"]
            for col in count_cols:
                intersect[f"{col}_max"] = intersect[col].apply(
                    lambda x: max([int(ac) for ac in x.split("; ")])
                )
                cols.append(f"{col}_max")
        # get max allele counts for TG
        except KeyError:
            count_cols = ["TG_AC", "seen_in_TG_count", "TG_nhomalt"]
            for col in count_cols:
                intersect[f"{col}_max"] = intersect[col].apply(
                    lambda x: max([int(ac) for ac in x.split("; ")])
                )
                cols.append(f"{col}_max")

    # merge population AF dataframe with annotSV df
    annotsv_pop_svs = pd.merge(
        annotsv_df,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={col: 0 for col in cols})
    return annotsv_pop_svs


def annotate_UCSC(chr, pos, end):
    UCSC_base_URL = (
        '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position='
    )
    UCSC_full_URL = f'{UCSC_base_URL}{chr}:{pos}-{end}","UCSC_link")'
    return UCSC_full_URL


def annotate_pb_regions(annotsv_df, regions, region_name):
    """
    Annotate SVs against PacBio odd regions or PacBio dark regions (bed files where fourth column indicates region affected)
    """
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    regions = pd.read_csv(regions, sep="\t")
    regions_bed = BedTool.from_dataframe(regions)
    intersect = annotsv_bed.intersect(
        regions_bed,
        wa=True,
        wb=True,
    ).to_dataframe()
    intersect.columns = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "ID",
        "CHROM_region",
        "POS_region",
        "END_region",
        region_name,
    ]
    # make a column with region details, e.g 1:25266309-25324509
    region_details_name = ("_").join(region_name.split("_")[:-1])
    intersect[region_details_name] = intersect[
        ["CHROM_region", "POS_region", "END_region"]
    ].apply(lambda x: f"{x[0]}:{x[1]}-{x[2]}", axis=1)
    # calculate percent of sample SV overlapped by region
    intersect[f"{region_details_name}_perc_overlap"] = intersect[
        ["POS", "END", "POS_region", "END_region"]
    ].apply(lambda x: calculate_sample_SV_overlap(x[0], x[1], x[2], x[3]), axis=1)
    intersect = intersect[
        [
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
            "ID",
            region_name,
            region_details_name,
            f"{region_details_name}_perc_overlap",
        ]
    ]
    cols = [col for col in intersect.columns if "PB" in col]
    # merge PB region dataframe with annotSV df
    intersect = intersect.astype(str)
    annotsv_odd_region_svs = pd.merge(
        annotsv_df,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={col: "." for col in cols})
    return annotsv_odd_region_svs


def annotate_repeats(annotsv_df, repeats):
    #  sample SV must be fully encompassed by repeat
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    repeats = pd.read_csv(repeats, sep="\t", header=None, names=["CHROM", "POS", "END", "Adotto_tandem_repeat"])
    repeats["CHROM"] = repeats["CHROM"].astype(str).apply(lambda x: x.replace("chr", ""))
    repeats_bed = BedTool.from_dataframe(repeats)
    intersect = annotsv_bed.intersect(
        repeats_bed,
        wa=True,
        wb=True,
        f=1,
    ).to_dataframe()
    intersect.columns = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "ID",
        "CHROM_repeat",
        "POS_repeat",
        "END_repeat",
        "Adotto_tandem_repeat",
    ]
    intersect = intersect[
        [
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
            "ID",
            "Adotto_tandem_repeat",
        ]
    ]

    # merge repeat dataframe with annotSV df
    intersect = intersect.astype(str)
    annotsv_repeat_svs = pd.merge(
        annotsv_df,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={"Adotto_tandem_repeat": "."})
    return annotsv_repeat_svs


def vcf_to_df(vcf_path):
    with open(vcf_path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    return pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})


def parse_snpeff(snpeff_df):
    svtype_list = []
    end_list = []
    ann_list = []
    pos_list = []
    # parse out svtype, end coordinate, ANN, and position
    for index, row in snpeff_df.iterrows():
        info = row["INFO"].split(";")
        pos = row["POS"]
        svtype = [i for i in info if "SVTYPE=" in i][0].split("=")[1]
        svtype_list.append(svtype)
        try:
            end = [i for i in info if "END=" in i][0].split("=")[1]
            if svtype == "BND" or svtype == "INS":
                end = int(end) + 1
        except IndexError:
            end = int(pos) + 1
        end_list.append(end)
        try:
            ann = [i for i in info if "ANN=" in i and "SVANN=" not in i][0]
            ann_list.append(ann)
        except IndexError:
            ann_list.append("NA")
        if svtype == "BND":
            CIPOS = [i for i in info if "CIPOS=" in i][0].split("=")[1].split(",")[0]
            pos = int(pos) + int(CIPOS)
        pos_list.append(pos)
    snpeff_df["SVTYPE"] = svtype_list
    snpeff_df["END"] = end_list
    snpeff_df["ANN"] = ann_list
    snpeff_df["POS"] = pos_list
    # from snpeff ANN field, parse out variant type, impact level (e.g. HIGH), gene name, and Ensembl gene ID
    variant_list_all = []
    impact_list_all = []
    gene_list_all = []
    ens_gene_list_all = []
    for line in snpeff_df["ANN"].values:
        variant_list = []
        impact_list = []
        gene_list = []
        ens_gene_list = []
        line = line.split(";")[-1]
        for anno in line.split(","):
            if "LOF=" in anno or "NMD=" in anno:
                pass
            else:
                try:
                    anno = anno.split("|")
                    variant_list.append(anno[1])
                    impact_list.append(anno[2])
                    # for some BNDs, snpeff doesn't annotate against a gene?
                    gene = anno[3]
                    if gene == "":
                        # add transcript instead
                        gene = anno[6]
                    gene_list.append(gene)
                    # same as above
                    ens_gene = anno[4]
                    if ens_gene == "":
                        ens_gene = anno[6]
                    ens_gene_list.append(ens_gene)
                except:
                    for l in [variant_list, impact_list, gene_list, ens_gene_list]:
                        l.append("NA")
        variant_list_all.append(";".join(set(variant_list)))
        impact_list_all.append(";".join(set(impact_list)))
        gene_list_all.append(";".join(set(gene_list)))
        ens_gene_list_all.append(";".join(set(ens_gene_list)))
    snpeff_df["IMPACT"] = impact_list_all
    snpeff_df["GENE_NAME"] = gene_list_all
    snpeff_df["ENSEMBL_GENE"] = ens_gene_list_all
    snpeff_df["VARIANT"] = variant_list_all
    snpeff_df["CHROM"] = [chr.replace("chr", "") for chr in snpeff_df["CHROM"].values]
    snpeff_df = snpeff_df[
        [
            "CHROM",
            "POS",
            "END",
            "ID",
            "SVTYPE",
            "ANN",
            "VARIANT",
            "IMPACT",
            "GENE_NAME",
            "ENSEMBL_GENE",
        ]
    ].astype(str)
    return snpeff_df


def merge_annotsv_snpeff(annotsv_df, snpeff_df):
    merged = annotsv_df.merge(
        snpeff_df, on=["CHROM", "POS", "END", "SVTYPE", "ID"], how="left"
    )
    merged = merged.drop(columns=["ANN"])
    return merged


def calculate_sample_SV_overlap(sample_pos, sample_end, database_pos, database_end):
    sample_len = sample_end - sample_pos
    overlap_start = max(sample_pos, database_pos)
    overlap_end = min(sample_end, database_end)
    overlap = (overlap_end - overlap_start) / float(sample_len) * 100
    overlap_perc = round(float(overlap), sigfigs=3)
    return overlap_perc


def add_clingen(clingen_df, gene, colname):
    genes = gene.split(";")
    clingen = []
    for gene in genes:
        # split by - for intergenic variants, which are annotated as <upstream_gene>-<downstream_gene>
        gene = gene.split("-")
        for g in gene:
            try:
                clingen.append(
                    str(clingen_df[clingen_df["Gene"] == g][colname].values[0])
                )

            except IndexError:
                pass
    if len(clingen) != 0:
        clingen = ",".join(clingen)
    else:
        clingen = "."
    return clingen


def add_BND_structure(svtype, info, alt):
    if svtype == "BND":
        info_extended = ";".join([info, alt])
        return info_extended
    else:
        return info


def main(
    df,
    snpeff_df,
    omim,
    hpo,
    vcf,
    prefix,
    exon_bed,
    gnomad,
    inhouse,
    tg,
    colorsdb,
    dark_regions,
    odd_regions,
    repeats,
    c4r,
    clingen_HI,
    clingen_TS,
    clingen_disease,
):
    print(c4r)
    # filter out SVs < 50bp
    df_len = df[apply_filter_length(df)]
    df_len = df_len.astype(str)
    # merge full and split AnnotSV annos
    df_merge = merge_full_split_annos(df_len)
    sample_cols = [col for col in df.columns if "PB" in col]
    if len(sample_cols) == 0:
        # for TCAG sequence IDs
        regexp = re.compile("[0-9][0-9]-[0-9]+")
        sample_cols = [col for col in df.columns if re.match(regexp, col)]
    if len(sample_cols) == 0:
        # C4R TCAG IDs
        sample_cols = [col for col in df.columns if "RLG" in col]
    if len(sample_cols) == 0:
        # genesteps TCAG IDs
        sample_cols = [col for col in df.columns if "RGS" in col]
    if len(sample_cols) == 0:
        # C4R IDs
        regexp = re.compile("\d+[A-Z]*_[A-Z]*\d+")
        sample_cols = [col for col in df.columns if re.match(regexp, col)]
    if len(sample_cols) == 0:
        # DECODER IDs
        sample_cols = [col for col in df.columns if "SK" in col]
    if len(sample_cols) == 0:
        # genoderm IDs
        for col in df.columns:
            print(col)
        sample_cols =  [col for col in df.columns if "GD" in col]
    if len(sample_cols) == 0:
       # ataxia IDS
       sample_cols =  [col for col in df.columns if "GYM" in col or "HSC" in col]
    if len(sample_cols) == 0:
        print("no sample cols identified")
        sys.exit(1)
                     
    # extract genotype and alt allele depth
    for sample in sample_cols:
        df_merge[f"{sample}_zyg"] = [
            get_genotype(row[sample])[0] for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_GT"] = [
            get_genotype(row[sample])[1] for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_AD"] = [
            get_alt_depth(row[sample]) for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_DP"] = [
            get_depth(row[sample]) for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_PS"] = [
            get_phase_set(row[sample]) for index, row in df_merge.iterrows()
        ]
    zyg_cols = [col for col in df_merge.columns if "_zyg" in col]
    gt_cols = [col for col in df_merge.columns if "_GT" in col]
    ad_cols = [col for col in df_merge.columns if "_AD" in col]
    dp_cols = [col for col in df_merge.columns if "_DP" in col]
    ps_cols = [col for col in df_merge.columns if "_PS" in col]

    # filter out benign SVs with AF > 1%
    # df_merge_notbenign = df_merge[apply_filter_benign(df_merge)]

    # add snpeff annos
    snpeff_df = parse_snpeff(snpeff_df)

    # merge snpeff and annotsv df
    df_merge = merge_annotsv_snpeff(df_merge, snpeff_df)

    # add HPO terms by gene matching
    try:
        df_merge["HPO"] = [
            add_hpo(hpo, gene) for gene in df_merge["ENSEMBL_GENE"].values
        ]
    except:
        print("No HPO terms")

    # add OMIM phenos and inheritance by gene matching
    print("Preparing OMIM data")
    omim_df = prepare_OMIM(f"{omim}/genemap2.txt")

    df_merge["omim_phenotype"] = [
        add_omim(omim_df, gene)[0] for gene in df_merge["ENSEMBL_GENE"].values
    ]
    df_merge["omim_inheritance"] = [
        add_omim(omim_df, gene)[1] for gene in df_merge["ENSEMBL_GENE"].values
    ]

    # add gnomAD SVs
    gnomad_cols = [
        "gnomad_NAME",
        "gnomad_GRPMAX_AF",
        "gnomad_AC",
        "gnomad_HOM",
    ]
    df_merge = annotate_pop_svs(df_merge, gnomad, gnomad_cols)

    # add C4R inhouse db SV counts
    inhouse_cols = [
        "C4R_ID",
        "C4R_REF",
        "C4R_ALT",
        "C4R_AC",
        "seen_in_C4R",
        "seen_in_C4R_count",
    ]

    df_merge = annotate_pop_svs(df_merge, inhouse, inhouse_cols)
    inhouse_cols = [col for col in inhouse_cols if col != "C4R_ID"]

    # add TG inhouse db SV counts
    tg_cols = [
        "TG_ID",
        "TG_REF",
        "TG_ALT",
        "TG_AC",
        "TG_nhomalt",
        "seen_in_TG",
        "seen_in_TG_count",
    ]

    df_merge  = annotate_pop_svs(df_merge, tg, tg_cols)
    tg_cols = [col for col in tg_cols if col != "TG_ID"]

    # add CoLoRSdb SVs
    colorsdb_cols = ["CoLoRSdb_SVLEN", "CoLoRSdb_AF", "CoLoRSdb_AC", "CoLoRSdb_AC_Hemi", "CoLoRSdb_nhomalt"]
    df_merge = annotate_pop_svs(df_merge, colorsdb, colorsdb_cols)
    df_merge = df_merge.drop(columns=["CoLoRSdb_SVLEN"])
    colorsdb_cols = [col for col in colorsdb_cols if col != "CoLoRSdb_SVLEN"]

    # add exon counts
    df_merge = get_exon_counts(df_merge, exon_bed)

    # add UCSC genome browser URL
    df_merge["UCSC_link"] = [
        annotate_UCSC(chrom, pos, end)
        for chrom, pos, end in zip(df_merge["CHROM"], df_merge["POS"], df_merge["END"])
    ]

    # add PacBio dark regions
    df_merge = annotate_pb_regions(df_merge, dark_regions, "PB_dark_region_gene")

    # add PacBio odd regions
    df_merge = annotate_pb_regions(df_merge, odd_regions, "PB_odd_region_type")

    # add PacBio repeats used in repeat expansion finding tool
    df_merge = annotate_repeats(df_merge, repeats)

    # add clingen haploinsufficiency and triplosensitivity scores
    df_merge["clingen_HI"] = [
        add_clingen(clingen_HI, gene, "Score") for gene in df_merge["GENE_NAME"].values
    ]
    df_merge["clingen_TS"] = [
        add_clingen(clingen_TS, gene, "Score") for gene in df_merge["GENE_NAME"].values
    ]
    df_merge["clingen_disease"] = [
        add_clingen(clingen_disease, gene, "Disease")
        for gene in df_merge["GENE_NAME"].values
    ]
    df_merge["clingen_classification"] = [
        add_clingen(clingen_disease, gene, "Classification")
        for gene in df_merge["GENE_NAME"].values
    ]

    # define columns to be included in report
    hpo_cols = ["HPO"] if isinstance(hpo, pd.DataFrame) else []

    # exclude C4R counts if not part of C4R study
    if c4r != "True":
        inhouse_cols = [col for col in inhouse_cols if col != "seen_in_C4R"]

    # add BND directionality information to INFO column
    df_merge["INFO_extended"] = df_merge.apply(
        lambda row: add_BND_structure(row["SVTYPE"], row["INFO"], row["ALT"]), axis=1
    )
    df_merge = df_merge.drop(["INFO"], axis=1)
    df_merge = df_merge.rename(columns={"INFO_extended": "INFO"})

    report_columns = (
        [
            "CHROM",
            "POS",
            "END",
            "SVLEN",
            "SVTYPE",
            "ID",
            "INFO",
            "FILTER",
            "GENE_NAME",
            "ENSEMBL_GENE",
            "VARIANT",
            "IMPACT",
            "UCSC_link",
            "omim_phenotype",
            "omim_inheritance",
            "DDD_mode",
            "DDD_consequence",
            "DDD_disease",
            "clingen_disease",
            "clingen_classification",
            "clingen_HI",
            "clingen_TS",
        ]
        + hpo_cols
        + zyg_cols
        + gt_cols
        + ad_cols
        + dp_cols
        + ps_cols
        + ["Tx", "Frameshift", "EXONS_SPANNED", "Nearest_SS_type", "Dist_nearest_SS"]
        + gnomad_cols
        + inhouse_cols
        + tg_cols
        + colorsdb_cols
        + [
            "ExAC_delZ",
            "ExAC_dupZ",
            "ExAC_cnvZ",
            "ExAC_synZ",
            "ExAC_misZ",
            "ExAC_pLI",
            "CytoBand",
            "RE_gene",
            "TAD_coordinate",
            "ENCODE_experiment",
            "Repeat_type_left",
            "Repeat_type_right",
            "SegDup_left",
            "SegDup_right",
            "PB_dark_region_gene",
            "PB_dark_region",
            "PB_dark_region_perc_overlap",
            "PB_odd_region_type",
            "PB_odd_region",
            "PB_odd_region_perc_overlap",
            "Adotto_tandem_repeat",
            "ENCODE_blacklist_characteristics_left",
            "ENCODE_blacklist_characteristics_right",
        ]
    )

    df_merge = df_merge[report_columns]
    df_merge = df_merge.drop(columns=["C4R_REF", "C4R_ALT", "TG_REF", "TG_ALT"])
    df_merge["GENE_NAME"] = df_merge["GENE_NAME"].replace("nan", ".")
    df_merge["omim_phenotype"].fillna("nan", inplace=True)
    df_merge["omim_inheritance"].fillna("nan", inplace=True)
    df_merge = df_merge.drop_duplicates()
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    df_merge.to_csv(f"{prefix}.pbsv.csv", index=False)
    df_merge.to_csv(f"{prefix}.pbsv.{today}.csv", index=False)

    # make a dictionary of variants with SVLEN >=50 and BNDs from vcf
    # validate that the only SVs missing from annotSV merged and SV >=50 bp are on non-canonical chr
    vcf_variant_dict = {}
    i = 0
    for record in vcf.fetch():
        i += 1
        try:
            length = abs(float(record.info["SVLEN"][0]))
            if length >= 50:
                vcf_variant_dict[record.id] = length
        except KeyError:
            if record.info["SVTYPE"] == "BND":
                vcf_variant_dict[record.id] = 0

    vcf_ID = vcf_variant_dict.keys()
    annotSV_ID = list(set(df_merge["ID"]))
    print("SVs in pbsv that are not present in AnnotSV tsv:")
    print(vcf_ID - annotSV_ID)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a structural variant report using an AnnotSV-annotated pbsv VCF generated from PacBio HiFi WGS"
    )
    parser.add_argument("-annotsv", type=str, help="AnnotSV tsv file", required=True)
    parser.add_argument("-snpeff", type=str, help="Snpeff vcf file", required=True)
    parser.add_argument(
        "-hpo",
        help="Tab delimited file containing gene names and HPO terms",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-omim",
        type=str,
        required=True,
        help="Path to directory containing OMIM mim2gene and morbidmap files",
    )
    parser.add_argument(
        "-vcf",
        help="pbsv vcf file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-exon",
        help="exon bed file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-gnomad",
        help="gnomad SVs in bed format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-inhouse",
        help="C4R inhouse database",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-tg",
        help="TG inhouse database",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-colorsdb",
        help="CoLoRSdb SVs in bed format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-dark_regions",
        help="PacBio dark regions",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-odd_regions",
        help="PacBio odd regions",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-repeats",
        help="Adotto repeats",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c4r",
        help="C4R sample?",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_HI",
        help="clingen HI scores",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_TS",
        help="clingen TS scores",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_disease",
        help="clingen disease association",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    # pull chrom, pos, end, SVtype, format fields and DDD annotations from AnnotSV text file
    df = pd.read_csv(args.annotsv, sep="\t", low_memory=False)
    # pull gene annotations from SnpEff VCF
    snpeff_df = vcf_to_df(args.snpeff)
    df = rename_SV_cols(df)
    prefix = args.annotsv.replace(".AnnotSV.tsv", "")
    # parse original pbsv VCF just to ensure all SVs are included in final report
    vcf = VariantFile(args.vcf)

    hpo = args.hpo
    if hpo:
        hpo = pd.read_csv(
            args.hpo,
            comment="#",
            skip_blank_lines=True,
            sep="\t",
            encoding="ISO-8859-1",
            engine="python",
        )
        # Phenotips TSV has a space in column name: " Gene symbol"
        hpo.columns = hpo.columns.str.strip()

    clingen_HI = pd.read_csv(
        args.clingen_HI,
        sep="\t",
        comment="#",
        header=None,
        names=["CHROM", "POS", "END", "Gene", "Score"],
    )
    clingen_TS = pd.read_csv(
        args.clingen_TS,
        sep="\t",
        comment="#",
        header=None,
        names=["CHROM", "POS", "END", "Gene", "Score"],
    )
    clingen_disease = pd.read_csv(args.clingen_disease, comment="#")

    main(
        df,
        snpeff_df,
        args.omim,
        hpo,
        vcf,
        prefix,
        args.exon,
        args.gnomad,
        args.inhouse,
        args.tg,
        args.colorsdb,
        args.dark_regions,
        args.odd_regions,
        args.repeats,
        args.c4r,
        clingen_HI,
        clingen_TS,
        clingen_disease,
    )
