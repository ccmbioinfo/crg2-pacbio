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
        self.df = pd.DataFrame(columns=columns, dtype=object)
        self._group_sv(all_sv)
        self.bedtool = self.make_ref_bedtool()

        for name in self.sample_list:
            self.df["%s_SV_DETAILS" % name] = list2string(self.df["%s_SV_DETAILS" % name])

        # append annotation fields to final df
        print(ann_df.head())
        print(self.df.head())
        ann_df = ann_df.set_index(keys=["CHROM", "START", "END", "SVTYPE"])
        self.df = self.df.join(ann_df, how='left')

        # retrieve genotype and CN fields (won't appear in joined dataframe for all samples if CNV coordinates are different between samples)


    def _parse_reports(self, report_paths):
        '''
            Merge all SV interval data from multiple vcf's in to a single BedTool instance

            Implementation:
                Use Panda's dataframe for some easy preprocessing, then create a BedTool from a tuple containing each row
        '''

        intervals = []
        sample_names = []
        ann_dfs = []

        sample_sv_fields =['CHROM', 'START', 'END', 'SVTYPE', 'samples', 'GT', 'CN']

        for report_path in report_paths:
            df = pd.read_csv(report_path, sep="\t") #use read_vcf because genotype field is not picked up with vcf_to_dataframe
            gt_col = [col for col in df.columns if '|GT' in col][0]
            cn_col = [col for col in df.columns if '|CN' in col][0]
            df.rename({gt_col: "GT", cn_col: "CN"}, axis=1, inplace=True)
            name = gt_col.split("|")[0]
            sample_names.append(name)
            df['samples'] = name

            intervals.extend(df[sample_sv_fields].itertuples(index=False))
            ann_dfs.append(df)

        ann_df = pd.concat(ann_dfs).astype(str)
        ann_df = ann_df.drop_duplicates(keep='first') #annotations for the same SV in a vcf can have slighly differing fields (ex. SVSCORE_MEAN)

        return BedTool(intervals), ann_df, sample_names

    def _group_sv(self, bedtool, reciprocal_overlap=0.5):
        already_grouped_intervals = set()

        for l in bedtool.intersect(bedtool, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap):

            ref_chr, ref_start, ref_end, ref_svtype, ref_name, ref_gt, ref_cn, \
            samp_chr, samp_start, samp_end, samp_svtype, samp_name, samp_gt, samp_cn = l

            ref_interval = (ref_chr, ref_start, ref_end, ref_svtype)
            samp_interval = (samp_chr, samp_start, samp_end, samp_svtype, samp_name, samp_gt, samp_cn)

            if (samp_interval not in already_grouped_intervals) and (ref_svtype == samp_svtype):
                self._add_interval(ref_interval, samp_interval)
                already_grouped_intervals.add(samp_interval)
        

    def _add_interval(self, ref_interval, samp_interval):
        samp_chr, samp_start, samp_end, samp_svtype, samp_name, samp_gt, samp_cn = samp_interval

        #Get reference to row
        ref_interval_short = ref_interval[0:4]
        try:
            self.df = self.df.set_index(keys=["CHROM", "START", "END", "SVTYPE"])
        except: 
            pass
        if ref_interval_short not in self.df.index:
            #make new row
            self.df.loc[ref_interval_short, :] = np.nan
            row = self.df.loc[ref_interval_short, :]
            row['N_SAMPLES'] = 0

            for name in self.sample_list:
                row[name] = 0
                row["%s_SV_DETAILS" % name] = []
        else:
            row = self.df.loc[ref_interval_short, :]

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

        print("adding interval")
        print(samp_gt, samp_cn)
        row["%s_SV_DETAILS" % samp_name].append('{}:{}-{}:{}:{}:{}'.format(samp_chr, samp_start, samp_end, samp_svtype, samp_gt, samp_cn))
    
    def write(self, outfile_name):
        self.df.to_csv(outfile_name, encoding='utf-8')

    def make_ref_bedtool(self):
        return BedTool(list(self.df.index.values))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merges CNV reports from TCAG to a single family report')
    parser.add_argument('-i', type=str, nargs='+', help='Report files containing CNV coordinates and annotations', required=True)
    parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.csv', required=True, type=str)
    args = parser.parse_args()
  

    cnvs = CNVGrouper(args.i)
    cnvs.df = cnvs.df.drop(['GT', 'CN'], axis=1)
    EXCEL_MAX = 32000
    for col in cnvs.df.columns:
        cnvs.df[col] = cnvs.df[col].apply(lambda x: str(x)[:EXCEL_MAX-1]) # truncate column contents
        cnvs.df[col] = ['.' if val == 'nan' else val for val in cnvs.df[col].tolist()]
    cnvs.write(args.o)
