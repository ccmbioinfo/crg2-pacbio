rule nanoplot:
    input:
        bam = get_bam
    output:
        "qc/nanoplot/{family}_{sample}/NanoStats.txt",
        "qc/nanoplot/{family}_{sample}/Non_Weighted_Histogram_Read_Length_mqc.png"
    log:
        "logs/nanoplot/{family}_{sample}.log"
    params:
        out_dir = "qc/nanoplot/{family}_{sample}"
    wrapper:
        get_wrapper_path("nanoplot")

rule qualimap:
    input: 
        bam = get_bam 
    output:
        "qc/qualimap/{family}_{sample}/genome_results.txt",
        "qc/qualimap/{family}_{sample}/qualimapReport.html",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt"
    params:
        out_dir = "qc/qualimap/{family}_{sample}",
        nw = config["params"]["qualimap"]["nw"],
        hm = config["params"]["qualimap"]["hm"],
        c = config["params"]["qualimap"]["c"],
        mem_size = config["params"]["qualimap"]["mem"],
        pipeline = config["run"]["pipeline"],
        extra = config["params"]["qualimap"]["extra"]
    threads: 8
    log:
        "logs/qualimap/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("qualimap")       

rule multiqc:
    input:
        [expand(input_file, sample=samples.index,family=project) for input_file in ["qc/nanoplot/{family}_{sample}/NanoStats.txt",
                                                                    "qc/qualimap/{family}_{sample}/genome_results.txt",                                                       
                                                                    "qc/qualimap/{family}_{sample}/qualimapReport.html",
                                                                    "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt",
                                                                    "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
                                                                    "qc/qualimap/{family}_{sample}/raw_data_qualimapReport/coverage_histogram.txt",                                 
                                                                    ]]                                                                                                                    
    output:
        report = "qc/multiqc/multiqc.html"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")
