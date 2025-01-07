rule define_trio_members:
    input: config["run"]["ped"] 
    output: "TRGT_denovo_IDs.txt"
    run:
        import pandas as pd
        # Load the pedigree file
        pedigree = pd.read_csv(input[0], sep=" ", header=None, names=["family_ID", "individual_ID", "paternal_ID", "maternal_ID", "sex", "phenotype"])
        print(pedigree)
        # Infer the roles of each sample
        children = pedigree[pedigree["paternal_ID"] != "0"][pedigree["maternal_ID"] != "0"]["individual_ID"].values
        with open(output[0], "w") as f:
            i = 0
            for child in children: # Could be multiple children. Phenotype field is not regularly populated, so just run pipeline on all children.
                father = pedigree[pedigree["individual_ID"] == pedigree[pedigree["individual_ID"] == child]["paternal_ID"].values[0]]["individual_ID"].values[0]
                mother = pedigree[pedigree["individual_ID"] == pedigree[pedigree["individual_ID"] == child]["maternal_ID"].values[0]]["individual_ID"].values[0]
                child = child.split("_")[1]
                mother = mother.split("_")[1]
                father = father.split("_")[1]
                if i == 0:
                    f.write(f"child\t{child}\n")
                    f.write(f"father\t{father}\n")
                    f.write(f"mother\t{mother}\n")
                    i += 1
                else:
                    f.write(f"child\t{child}\n")

rule trgt_denovo:
    input: 
        ID_map = "TRGT_denovo_IDs.txt",
        samples = config["run"]["samples"]
    params:
        trgt_denovo = config["tools"]["trgt-denovo"],
        ref = config["ref"]["genome"],
        bed = config["trgt"]["adotto_repeats"]
    output: "TRGT_denovo/{family}_{child}.TRGT.denovo.tsv"
    log: "logs/denovo_TRs/{family}_{child}.TRGT-denovo.log"
    resources:
        threads = 8
    shell:
        """
        # Get family member IDs from mapping file
        while IFS=$'\t' read -r family_member sample_ID
        do
            if [[ "${{family_member}}" == "father" ]]; then
                father_ID=${{sample_ID}}
            elif [[ "${{family_member}}" == "mother" ]]; then
                mother_ID=${{sample_ID}}
            fi
        done < {input.ID_map}

        # Get BAM paths
        child_bam=`cat {input.samples} | grep -P "{wildcards.child}\t" | awk '{{print $2}}' | sed 's/.bam/.trgt/'`
        father_bam=`cat {input.samples} | grep -P "${{father_ID}}\t" | awk '{{print $2}}' | sed 's/.bam/.trgt/'`
        mother_bam=`cat {input.samples} | grep -P "${{mother_ID}}\t" | awk '{{print $2}}' | sed 's/.bam/.trgt/'`
	echo         {params.trgt_denovo} trio --reference {params.ref} \
                --bed {params.bed} \
                --father $father_bam \
                --mother $mother_bam \
                --child $child_bam \
                -@ {resources.threads} \
                --out {output}

        {params.trgt_denovo} trio --reference {params.ref} \
                --bed {params.bed} \
                --father $father_bam \
                --mother $mother_bam \
                --child $child_bam \
                -@ {resources.threads} \
                --out {output}
        """   

rule annotate_trgt_denovo:
    input: "TRGT_denovo/{family}_{child}.TRGT.denovo.tsv"
    output: "TRGT_denovo/{family}_{child}.TRGT.denovo.annotated.csv"
    params:
      crg2_pacbio = config["tools"]["crg2_pacbio"],
      genes = config["trgt"]["ensembl"],
      constraint = config["trgt"]["gnomad_constraint"],
      OMIM = config["annotation"]["omim_path"],
      segdup = config["trgt"]["segdup"],
      controls = config["trgt"]["control_alleles"],
      c4r = config["annotation"]["c4r"],
      HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
      c4r_outliers = config["trgt"]["C4R_outliers"]
    log:  "logs/denovo_TRs/{family}_{child}.annotate.TRGT.denovo.log"
    conda: 
        "../envs/str_sv.yaml"
    shell:
        """
        if [[ {params.HPO} == "none" ]]
        then
            (python3 {params.crg2_pacbio}/scripts/annotate_denovo_TRs.py --repeats {input} \
                --output_file  {output} \
                --ensembl {params.genes} \
                --gnomad_constraint {params.constraint} \
                --OMIM_path {params.OMIM} \
                --segdup {params.segdup} \
                --c4r {params.c4r} \
                --controls {params.controls} \
		        --c4r_outliers {params.c4r_outliers}) > {log} 2>&1
        else
            (python3 {params.crg2_pacbio}/scripts/annotate_denovo_TRs.py --repeats {input} \
                --output_file  {output} \
                --ensembl {params.genes} \
                --gnomad_constraint {params.constraint} \
                --OMIM_path {params.OMIM} \
                --segdup {params.segdup} \
                --c4r {params.c4r} \
                --controls {params.controls} \
                --hpo {params.HPO} \
		        --c4r_outliers {params.c4r_outliers}) > {log} 2>&1
        fi
        """
