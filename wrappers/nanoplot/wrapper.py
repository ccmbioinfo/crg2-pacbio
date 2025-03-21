from snakemake.shell import shell

out_dir = snakemake.params.get("out_dir") 

shell(
    "NanoPlot -t 4 --bam {snakemake.input} --N50 --plots hex --outdir {out_dir}"
    "&& mv {out_dir}/Non_weightedHistogramReadlength.png {out_dir}/Non_Weighted_Histogram_Read_Length_mqc.png"
    )
