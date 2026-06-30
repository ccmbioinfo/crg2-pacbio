from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
tool = snakemake.params.tool
crg2_pacbio = snakemake.params.crg2_pacbio
input_vcf = snakemake.input[0]
pythonpath = tool.replace("bin", "")
python = " export PYTHONPATH={pythonpath}; "

# Use the bundled mity chrM reference to match mity's normalization/report expectations.
reference_fasta = pythonpath + "/mitylib/reference/hg38.chrM.fa"

bcftools_norm = (
    " bcftools norm -f {reference_fasta} -m-both -Oz "
    "-o {outdir}/{family}.mt.normalise.vcf.gz {input_vcf};"
)
vt = (
    " vt decompose_blocksub -o {outdir}/{family}.mt.normalise.decompose.raw.vcf.gz "
    "{outdir}/{family}.mt.normalise.vcf.gz;"
)
fix_vcf = (
    " python {crg2_pacbio}/scripts/format_mitorsaw_mt_vcf.py "
    "{outdir}/{family}.mt.normalise.decompose.raw.vcf.gz "
    "{outdir}/{family}.mt.normalise.decompose.vcf.gz;"
)
tabix = " tabix -f {outdir}/{family}.mt.normalise.decompose.vcf.gz;"
remove_intermediate = (
    " rm {outdir}/{family}.mt.normalise.vcf.gz "
    "{outdir}/{family}.mt.normalise.decompose.raw.vcf.gz"
)

shell("(" + python + bcftools_norm + vt + fix_vcf + tabix + remove_intermediate + ") {log}")
