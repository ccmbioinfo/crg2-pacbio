#!/usr/bin/env python3
"""
Fix CIPOS field for CNV variants at chromosome start (position 1).

For variants with position 1, modify the CIPOS field to ensure it doesn't
have negative values, since the position cannot be negative.
"""

import pysam
import sys
import argparse


def fix_cipos_for_chr_start(input_vcf, output_vcf):
    """
    Read VCF file and modify CIPOS for variants at position 1.
    
    Args:
        input_vcf: Path to input VCF file
        output_vcf: Path to output VCF file
    """
    with pysam.VariantFile(input_vcf, 'r') as vcf_in:
        # Create output VCF with same header
        with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
            
            for variant in vcf_in:
                # Check if variant is at position 1 (chromosome start)
                if variant.pos == 1 and 'CIPOS' in variant.info:
                    # Get current CIPOS value
                    cipos = variant.info['CIPOS']
                    
                    # CIPOS can be a tuple (min, max) or a list
                    if isinstance(cipos, (tuple, list)) and len(cipos) == 2:
                        # Ensure minimum value is not negative
                        cipos_min = max(0, cipos[0])
                        cipos_max = cipos[1]
                        variant.info['CIPOS'] = (cipos_min, cipos_max)
                    elif isinstance(cipos, (tuple, list)) and len(cipos) == 1:
                        # Single value case
                        cipos_val = max(0, cipos[0])
                        variant.info['CIPOS'] = (cipos_val,)
                        modified_count += 1
                    elif isinstance(cipos, (int, float)):
                        # Single numeric value
                        variant.info['CIPOS'] = (max(0, int(cipos)),)
                        modified_count += 1
                
                # Write variant (modified or unmodified) to output
                vcf_out.write(variant)
    

def main():
    parser = argparse.ArgumentParser(
        description='Fix CIPOS field for CNV variants at chromosome start (position 1)'
    )
    parser.add_argument(
        '-input_vcf',
        help='Input VCF file'
    )
    parser.add_argument(
        '-output_vcf',
        help='Output VCF file'
    )
    
    args = parser.parse_args()
    
    fix_cipos_for_chr_start(args.input_vcf, args.output_vcf)


if __name__ == '__main__':
    main()

