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
    output: "TRGT_denovo/{family}_TRGT_denovo.tsv"
    resources:
        threads = 8
    shell:
        """
        # map family member label (e.g. child) to sample ID
        while IFS=$'\t' read -r family_member sample_ID
        do
            echo $sample_ID
            if [[ "${{family_member}}" == "child" ]]; then
                child_ID=`echo ${{sample_ID}}`
            elif [[ "${{family_member}}" == "father" ]]; then
                father_ID=`echo ${{sample_ID}}`
            else
                mother_ID=`echo ${{sample_ID}}`
            fi
        done<{input.ID_map}

        child_bam=`awk '{{print $2}}' {input.samples} | grep ${{child_ID}} | sed 's/.bam/.trgt/'`
        father_bam=`awk '{{print $2}}' {input.samples} | grep ${{father_ID}} | sed 's/.bam/.trgt/'`
        mother_bam=`awk '{{print $2}}' {input.samples} | grep ${{mother_ID}} | sed 's/.bam/.trgt/'`

        {params.trgt_denovo} trio --reference {params.ref} \
                --bed {params.bed} \
                --father $father_bam \
                --mother $mother_bam \
                --child $child_bam \
                -@ {resources.threads} \
                --out {output}
                                                            
        """    

