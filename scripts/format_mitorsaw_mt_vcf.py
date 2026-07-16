import sys

import pysam


input_vcf = sys.argv[1]
output_vcf = sys.argv[2]


original_vcf = pysam.VariantFile(input_vcf)
formatted_vcf = pysam.VariantFile(output_vcf, "wz", header=original_vcf.header)

for record in original_vcf:
    if not list(record.filter.keys()):
        record.filter.add("PASS")

    for sample in record.samples.values():
        if "VAF" not in record.format.keys():
            continue

        vaf = sample.get("VAF", None)
        if vaf is None:
            sample["VAF"] = (0.0,)
        else:
            sample["VAF"] = tuple(0.0 if x is None else x for x in vaf)

    formatted_vcf.write(record)

original_vcf.close()
formatted_vcf.close()
