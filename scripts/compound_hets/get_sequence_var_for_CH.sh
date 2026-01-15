#!/bin/bash
#  exports gemini.db database to gemini.db.txt file
#  database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table


if [ -z $file ]
then
    file=$1
fi

severity_threshold=$2
max_af=0.01
alt_depth=3

gemini query -q "select name from samples order by name" $file > samples.txt


sQuery="select \
        chrom as Chrom,\
        start+1 as Pos,\
        variant_id as Variant_id,\
        ref as Ref,\
        alt as Alt,\
        impact as Variation,\
        dp as Depth,\
        qual as Quality,\
        gene as Gene,\
        COALESCE(clinvar_pathogenic, '') || COALESCE( ';' || NULLIF(clinvar_sig,''), '') || COALESCE( ';' || NULLIF(clinvar_sig_conf,''), '') as Clinvar, \
        ensembl_gene_id as Ensembl_gene_id,\
        gnomad_af_grpmax as Gnomad_af_grpmax,\
        tg_lrwgs_ac as TG_LRWGS_ac,\
        cadd_phred as Cadd_score,\
        COALESCE(spliceai_score, '') as SpliceAI_score,
        promoterAI as promoterAI_score,
        PS as PS,"

while read sample
do
	sQuery=$sQuery"gts."$sample","
    sQuery=$sQuery"gt_types."$sample","
    sQuery=$sQuery"gt_phases."$sample","
	sQuery=$sQuery"gt_alt_depths."$sample","
	sQuery=$sQuery"gt_depths."$sample","
    sQuery=$sQuery"gt_quals."$sample","
done < samples.txt


if [[ "$severity_threshold" == 'HIGH-MED' ]]
then
    severity_filter=" and impact_severity<>'LOW' " # exclude LOW impact variants
else
    severity_filter=" and impact_severity == 'LOW' "
fi

sQuery=$sQuery"hgvsc as Nucleotide_change_ensembl,\
        hgvsp as Protein_change_ensembl from variants"

initialQuery=$sQuery # keep the field selection part for later use

sQuery=$sQuery" where gnomad_af_grpmax <= "${max_af}" "$caller_filter""${severity_filter}""


# keep variant where the alt depth is >=3 in any one of the samples or they're all -1 (sometimes happens for freebayes called variants?)
s_gt_filter="(gt_alt_depths).(*).(>="${alt_depth}").(any) or (gt_alt_depths).(*).(==-1).(all)"
gemini query -q "$sQuery" --gt-filter "${s_gt_filter}" --header $file

# also get the clinvar variants (duplicates will be removed later)
cQuery=$initialQuery
cQuery=$cQuery" where gnomad_af_grpmax <= ${max_af} "$caller_filter" and Clinvar <> ''"
# only get variants where AD >= 1 (any sample with an alternate read)
c_gt_filter="(gt_alt_depths).(*).(>=1).(any) or (gt_alt_depths).(*).(==-1).(all)"
gemini query -q "$cQuery" --gt-filter "$c_gt_filter" $file

# if allele frequency is > 1% and Clinvar is pathogenic, likely pathogenic or conflicting and any status except for no assertion
cQuery=$initialQuery
cQuery=$cQuery" where gnomad_af_grpmax > ${max_af} "$caller_filter" and Clinvar_status != 'no_assertion_criteria_provided' and Clinvar in ('Pathogenic', 'Likely_pathogenic', 'Conflicting_interpretations_of_pathogenicity')"
gemini query -q "$cQuery" $file
