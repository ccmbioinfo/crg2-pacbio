from snakemake.shell import shell

out_dir = snakemake.params.get("out_dir") 

shell("NanoPlot -t 4 --bam {snakemake.input} --N50 --no_static --plots hex --outdir {out_dir}")