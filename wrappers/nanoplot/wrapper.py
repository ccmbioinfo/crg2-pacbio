from snakemake.shell import shell

out_dir = snakemake.params.get("out_dir")
prefix = snakemake.params.get("prefix")

shell:

"NanoPlot -t 4 --bam {snakemake.input} --prefix {prefix} --N50 --no_static --plots hex --outdir {out_dir}
