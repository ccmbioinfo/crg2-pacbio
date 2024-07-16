import numpy as np
import pandas as pd
import os
import argparse
from pybedtools import BedTool

# source: https://github.com/ccmbioinfo/crg/blob/master/merge.cnv.reports.py

class CNVGrouper:
    def __init__(self, reports):
        def list2string(col):
            return col.apply(lambda x: ', '.join(x) if isinstance(x, list) else x)

        #key - dataframe col name : value - vcf field name
        self.index_cols = ['CHROM', 'START', 'END', 'SVTYPE']

        all_sv, ann_df, sample_list = self._parse_reports(reports)
        assert len(sample_list) == len(set(sample_list)), "Duplicate sample names among input vcf's detected: %s" % sample_list

        columns = ['CHROM', 'START', 'END', 'SVTYPE', 'N_SAMPLES']
        columns.extend(sample_list)
        columns.extend(["%s_SV_DETAILS" % s for s in sample_list])

        self.sample_list = sample_list
        self.df = pd.DataFrame(columns=columns).set_index(keys=self.index_cols)
        self._group_sv(all_sv)
        self.bedtool = self.make_ref_bedtool()

        for name in self.sample_list:
            self.df["%s_SV_DETAILS" % name] = list2string(self.df["%s_SV_DETAILS" % name])

        # append annotation fields to final df
        self.df = self.df.join(ann_df, how='left')

    def _parse_reports(self, report_paths):
        '''
            Merge all SV interval data from multiple vcf's in to a single BedTool instance

            Implementation:
                Use Panda's dataframe for some easy preprocessing, then create a BedTool from a tuple containing each row
        '''

        intervals = []
        sample_names = []
        ann_dfs = []

        sample_sv_fields = self.index_cols + ['samples']

        for report_path in report_paths:
            df = pd.read_csv(report_path, sep="\t") #use read_vcf because genotype field is not picked up with vcf_to_dataframe

            df = df.astype({'sample': 'string'})
            name = df.pop('sample')[0]
            sample_names.append(name)
            df['samples'] = name

            intervals.extend(df[sample_sv_fields].itertuples(index=False))
            ann_dfs.append(df)

        ann_df = pd.concat(ann_dfs).astype(str).set_index(self.index_cols)
        ann_df = ann_df[~ann_df.index.duplicated(keep='first')] #annotations for the same SV in a vcf can have slighly differing fields (ex. SVSCORE_MEAN)

        return BedTool(intervals), ann_df, sample_names

    def _group_sv(self, bedtool, reciprocal_overlap=0.5):
        already_grouped_intervals = set()

        for l in bedtool.intersect(bedtool, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap):

            ref_chr, ref_start, ref_end, ref_svtype, ref_name, \
            samp_chr, samp_start, samp_end, samp_svtype, samp_name = l

            ref_interval = (ref_chr, ref_start, ref_end, ref_svtype)
            samp_interval = (samp_chr, samp_start, samp_end, samp_svtype, samp_name)

            if (samp_interval not in already_grouped_intervals) and (ref_svtype == samp_svtype):
                self._add_interval(ref_interval, samp_interval)
                already_grouped_intervals.add(samp_interval)
        
        self.df.sort_index(inplace=True)

    def _add_interval(self, ref_interval, samp_interval):
        samp_chr, samp_start, samp_end, samp_svtype, samp_name = samp_interval

        #Get reference to row
        if ref_interval not in self.df.index:
            #make new row
            self.df.loc[ref_interval, :] = np.nan
            row = self.df.loc[ref_interval, :]
            row['N_SAMPLES'] = 0

            for name in self.sample_list:
                row[name] = 0
                row["%s_SV_DETAILS" % name] = []
        else:
            row = self.df.loc[ref_interval, :]

        #Set values for row
        try:
            if row[samp_name] == 0:
                row['N_SAMPLES'] += 1
                row[samp_name] = 1
        except:
            samp_name = int(samp_name)
            if row[samp_name] == 0:
                row['N_SAMPLES'] += 1
                row[samp_name] = 1

        row["%s_SV_DETAILS" % samp_name].append('{}:{}-{}:{}'.format(samp_chr, samp_start, samp_end, samp_svtype))
    
    def write(self, outfile_name):
        self.df.to_csv(outfile_name, sep='\t', encoding='utf-8')

    def make_ref_bedtool(self):
        return BedTool(list(self.df.index.values))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merges CNV reports from TCAG to a single family report')
    parser.add_argument('-i', type=str, nargs='+', help='Report files containing CNV coordinates and annotations', required=True)
    parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.tsv', required=True, type=str)
    args = parser.parse_args()

    drop_cols = ["SIZE", "FILTER", "ALT", "Num_CNVs", "Length_CNVs", "Length_Gaps", "Percent_Gap", \
    "cnvn_count", "cnvn_details", "chrAnn", "startAnn", "endAnn", "variantTypeAnn", "CNVsizeAnn", "ASDgenes", \
    "ISCA_region", "CNV_ISCA_percOverlap", "ISCA_CNV_percOverlap", "ISCA_CNV_percOverlap_max", \
    "ISCA_matchCNVmax_percOverlap", "exon_symbol_ISCA", "DGVpercFreq_subjects_coverageStudies_50percRecipOverlap", \
    "fileID", "comment", "samples", "CNVnumberOfGeneSymbols", \
    "gnomAD_commonSV", "gnomAD_rareSV", "gnomAD_oe_lof", "gnomAD_oe_lof_upper", "gnomAD_oe_mis", "gnomAD_oe_mis_upper", "gnomAD_pLI", "gnomAD_pRec", ]

    end_cols = ["SAMPLE|REFCN", "SAMPLE|CN", "IMPRECISE", "cnv_type_conflict", "cnv_type_confict_coverage",\
    "MPO_NervousSystem", "MPO_Growth", "MPO_Other", "HPO_NervousSystem", "HPO_Growth", "HPO_Other", ]

    first_cols = ["gene_symbol", "gene_egID", "gene_symbol_CNVstart", "gene_symbol_CNVend", "exon_symbol", "exon_egID", \
    "cds_symbol", "cds_egID", "ExAC_pLI", "repeatMasker_percOverlap", "dirtyRegion_percOverlap", "CGD", "OMIM_MorbidMap", \
    "cnvnatorPercFreq_50percRecipOverlap", "erdsPercFreq_50percRecipOverlap", ]    

    cnvs = CNVGrouper(args.i)
    cnvs.df = cnvs.df.drop(columns=drop_cols)
    
    report_cols = [col for col in cnvs.df.columns if col not in set(end_cols + first_cols)]
    first_index = int([i for i, col in enumerate(report_cols) if '_SV_DETAILS' in col][-1]) + 1
    ordered_cols = report_cols[:first_index] + first_cols + report_cols[first_index:] + end_cols

    cnvs.df = cnvs.df[ordered_cols]
    cnvs.write(args.o)
