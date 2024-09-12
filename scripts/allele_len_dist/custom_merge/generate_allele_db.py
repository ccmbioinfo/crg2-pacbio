import argparse
import os
import pandas as pd
from pathlib import Path, PosixPath
import gzip

# https://github.com/egor-dolzhenko/discover-expansions/blob/main/discover-expansions.ipynb
# MC: added motif purity and methylation

def skip_header(file_handle: gzip.GzipFile, prefix: bytes = b'#') -> None:
    last_pos = file_handle.tell()
    while file_handle.readline().startswith(prefix):
        last_pos = file_handle.tell()
    file_handle.seek(last_pos)

def get_alleles(path: PosixPath):    
    gt_actions = {
        "0/0": lambda: (trid, motif_purity, avg_methylation, [ref, ref]),
        "0/1": lambda: (trid, motif_purity, avg_methylation, [ref, alt]),
        "1/2": lambda: (trid, motif_purity, avg_methylation, alt.split(",")),
        "1/1": lambda: (trid, motif_purity, avg_methylation, [alt, alt]),
        "1": lambda: (trid, motif_purity, avg_methylation, [alt]),
        "0": lambda: (trid, motif_purity, avg_methylation, [ref])
    }
    
    with gzip.open(path, 'r') as f_in:
        skip_header(f_in)
        for line in f_in:     
            line = line.decode("utf8")
            sl = line.split()
            
            gt = sl[-1].split(":")[0]
            
            if gt == '.':
                continue
            
            assert gt in gt_actions, f"Unknown gt: {gt}"
            
            ref, alt = sl[3], sl[4]               
            trid = sl[-3].split(";")[0].lstrip("TRID=")
            motif = sl[-3].split(";")[2].lstrip("MOTIFS=")
            motif_purity = sl[-1].split(":")[-2]
            avg_methylation = sl[-1].split(":")[-1]

            yield gt_actions[gt]()          

def main(vcf_path, out_filepath):     
    print(vcf_path)
    vcf_path = Path(vcf_path).resolve(strict=True) 
    with gzip.open(out_filepath, "w", compresslevel=6) as f_out:
        for index, path in enumerate(vcf_path.glob("*.vcf.gz")):    
            sample = path.name.rstrip(".vcf.gz")
            for (trid, motif_purity, avg_methylation, alleles) in get_alleles(path):        
                alleles = ",".join(alleles)
                f_out.write(f"{trid} {sample} {motif_purity} {avg_methylation} {alleles}\n".encode())  


if __name__ == "__main__":
    # if running from the command-line
    description = "Generate repeat allele database from TRGT VCF"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--vcf_path",
        type=str,
        required=True,
        help="Path to TRGT VCFs",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )

    args = parser.parse_args()
    print("Generating repeat allele database")
    main(args.vcf_path, args.output_file)