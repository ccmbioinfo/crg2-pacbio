import argparse
import os
import pandas as pd
from pathlib import Path, PosixPath
import gzip

# https://github.com/egor-dolzhenko/discover-expansions/blob/main/discover-expansions.ipynb

def skip_header(file_handle: gzip.GzipFile, prefix: bytes = b'#') -> None:
    last_pos = file_handle.tell()
    while file_handle.readline().startswith(prefix):
        last_pos = file_handle.tell()
    file_handle.seek(last_pos)

def get_alleles(path: PosixPath):    
    gt_actions = {
        "0/0": lambda: (trid, [ref, ref]),
        "0/1": lambda: (trid, [ref, alt]),
        "1/2": lambda: (trid, alt.split(",")),
        "1/1": lambda: (trid, [alt, alt]),
        "1": lambda: (trid, [alt]),
        "0": lambda: (trid, [ref])
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

            yield gt_actions[gt]()        

def main(vcf_path, out_filepath):     
    if not os.path.exists("repeat_outliers"):
        Path("repeat_outliers").mkdir(exist_ok=True)

    print(vcf_path)
    vcf_path = Path(vcf_path).resolve(strict=True) 
    with gzip.open(out_filepath, "w", compresslevel=6) as f_out:
        for index, path in enumerate(vcf_path.glob("*.vcf.gz")):    
            sample = path.name.rstrip(".vcf.gz")
            for (trid, alleles) in get_alleles(path):        
                alleles = ",".join(alleles)
                f_out.write(f"{trid} {sample} {alleles}\n".encode())  


if __name__ == "__main__":
    # if running from the command-line
    description = "Generate repeat allele database from TRGT VCF"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--vcf_path",
        type=str,
        required=True,
        help="Path to TRGT VCFs for cases and controls",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output filepath",
    )
    parser.add_argument(
        "--family",
        type=str,
        required=True,
        help="Family ID",
    )

    args = parser.parse_args()
    units = pd.read_table(args.vcf_path, dtype=str).set_index(["family"], drop=False)
    input_vcf_path=units.loc[args.family, "trgt_vcf_dir"]
    print("Generating repeat allele database")
    main(args.input_vcf_path, args.output_file)