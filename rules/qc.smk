rule nanoplot:
    input:
        bam = "mapped/{family}_{sample}.bam"
    output:
        "qc/nanoplot/{family}_{sample}/NanoStats.txt",
        "qc/nanoplot/{family}_{sample}/Non_Weighted_Histogram_Read_Length_mqc.png"
    log:
        "logs/nanoplot/{family}_{sample}.log"
    params:
        out_dir = "qc/nanoplot/{family}_{sample}"
    wrapper:
        get_wrapper_path("nanoplot")

rule multiqc:
    input:
        [expand(input_file, sample=samples.index,family=project) for input_file in ["qc/nanoplot/{family}_{sample}/NanoStats.txt" 
                                                                    ]]                                                                                                                    
    output:
        report = "qc/multiqc/multiqc.html"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        get_wrapper_path("multiqc")

