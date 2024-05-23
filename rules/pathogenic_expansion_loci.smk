rule genotype_pathogenic_loci:
    input: get_bam
    output: 
        vcf = temp("pathogenic_repeats/{family}_{sample}.trgt.unsorted.vcf.gz"),
        spanning_bam = "pathogenic_repeats/{family}_{sample}.trgt.unsorted.spanning.bam"
    params: 
        trgt = config["tools"]["trgt"],
        ref = config["ref"]["genome"],
        path_repeats = config["annotation"]["pathogenic_repeats"]["trgt_catalog"]
    log: "logs/pathogenic_repeats/{family}_{sample}.trgt.log"
    conda:
        "../envs/common.yaml"
    shell: 
        """
        # get sex, code from https://github.com/ccmbioinfo/crg/get_XY.sh
        x=`samtools idxstats {input} | egrep  "X|chrX"`;
        y=`samtools idxstats {input} | egrep  "Y|chrY"`;
        xcov=`echo $x | awk '{{ printf("%0.5f", $3/$2); }}'`;
        ycov=`echo $y | awk '{{ printf("%0.5f", $3/$2); }}'`;

        rat=$(echo "scale=4; ${{xcov}}/${{ycov}}" | bc)
        if (( $(echo "$rat > 5.0" | bc -l) )); then
            sex=XX
        else
            sex=XY
        fi

        echo "$sex"

        {params.trgt} genotype --genome {params.ref} \
            --reads {input} \
            --repeats {params.path_repeats} \
            --output-prefix pathogenic_repeats/{wildcards.family}_{wildcards.sample}.trgt.unsorted \
            --karyotype $sex
        """

rule sort_trgt_vcf:
    input: "pathogenic_repeats/{family}_{sample}.trgt.unsorted.vcf.gz"
    output: "pathogenic_repeats/{family}_{sample}.trgt.vcf.gz"
    log: "logs/bcftools/{family}_{sample}.sort.trgt.log"
    conda:
        "../envs/common.yaml"
    shell:
        """
        bcftools sort -Ob -o {output} {input};
        bcftools index {output}
        """

rule merge_vcfs:
    input: 
       vcfs=expand("pathogenic_repeats/{{family}}_{sample}.trgt.vcf.gz", sample=samples.index)
    output: "pathogenic_repeats/{family}.known.path.str.loci.vcf"
    log: "logs/pathogenic_repeats/{family}.merge.repeat.vcfs.log"
    conda:
        "../envs/common.yaml"
    shell: 
        """
        # determine if there is one TRGT VCF or multiple
        vcfs={input.vcfs} 
        count=0
        for file in ${{vcfs[@]}} ; do
        if [[ -f "$file" ]]; then
            count=$((count+1))
        fi
        done

        if [[ $count -eq 1 ]]; then
            echo "Single VCF"
            input_prefix=`echo {input.vcfs} | sed 's/.gz//'`
            gunzip {input.vcfs}
            mv $input_prefix {output}
        else
            echo "Multiple VCFs, merging"
            bcftools merge {input.vcfs} -o {output} -O v -F x -m all
        fi
        """

rule annotate_pathogenic_repeats:
    input: "pathogenic_repeats/{family}.known.path.str.loci.vcf"
    output: "pathogenic_repeats/{family}.known.path.str.loci.csv"
    params:
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        disease_thresholds = config["annotation"]["pathogenic_repeats"]["disease_thresholds"]
    log: "logs/pathogenic_repeats/{family}.annotate.repeats.log"
    resources:
        mem_mb = 10000
    conda:
        "../envs/str_sv.yaml"
    shell: 
        """
        (python3 {params.crg2_pacbio}/scripts/annotate_path_str_loci.py --vcf {input} \
            --output_file  {output} \
            --disease_thresholds {params.disease_thresholds}) > {log} 2>&1
        """