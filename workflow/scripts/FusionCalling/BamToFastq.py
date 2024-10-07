#!/usr/bin/env python3

# Adapted from https://github.com/TrinityCTAT/CTAT-LR-fusion/blob/main/util/sc/10x_ubam_to_fastq.py

import argparse
import pysam
import re

def bam_to_fastq(bam,fastq):

    with open(fastq, 'w') as f:
        samreader = pysam.AlignmentFile(bam, "rb", check_sq=False)
        for read in samreader:

            d = read.to_dict()
            
            read_name = d['name']
            read_seq = d['seq']
            quals = d['qual']

            # 10x tag descriptions at:
            # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam#:~:text=Barcoded%20BAM%20Tags,-The%20cellranger%20pipeline%20outputs%20an

            CB = "NA"
            if read.has_tag("CB"):
                cell_barcode = read.get_tag("CB", "Z")[0]
                cell_barcode = re.sub("-1$", "", cell_barcode)

            umi = "NA"
            if read.has_tag("UB"):
                umi = read.get_tag("UB", "Z")[0]
            elif len(read_name.split('.'))>=2:
                umi = read_name.split('.')[-2][:-3]
            

            
            read_name = "^".join([cell_barcode, umi, read_name])
            
            f.write("@" + read_name + '\n')
            f.write(read_seq + '\n')
            f.write("+" + '\n')
            f.write(quals + '\n')

def initialize_parser():
    parser = argparse.ArgumentParser(description='Rename cell types in cancer/non-cancer')
    parser.add_argument('--bam', type=str, default=1, help='User input barcode file', required = True)
    parser.add_argument('--fastq', type=str, default=1, help='Barcode file with redefined cancer/non-cancer celltypes', required = True)
    return (parser)


def main():
    # 1. Arguments
    parser = initialize_parser()
    args = parser.parse_args()

    bam=args.bam
    fastq=args.fastq

    # 2. Bam to fastq
    bam_to_fastq(bam,fastq)

#-------------
# Execute code
#-------------

if __name__ == '__main__':
    main()
