
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


fill_tags_cmd = (
    "bcftools +fill-tags {snakemake.input} " 
    "-o {snakemake.output} -O v "
    "-- -t 'DP:1=int(sum(FORMAT/DP))'"
)

shell(fill_tags_cmd + " {log}")
