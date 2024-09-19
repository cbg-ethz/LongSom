#!/cluster/work/bewi/members/dondia/Anaconda3/envs/pysam/bin/python

import pandas as pd
import numpy as np
import pysam
from pathlib import Path
import glob
import argparse
import gzip
from collections import defaultdict

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq)
    letters = [complement[base] for base in bases] 
    letters = ''.join(letters)
    reverse = letters[::-1]
    return reverse

def get_emptydrops(emptydroplets):
    emptydrops = pd.read_csv(emptydroplets, sep = '\t').barcodes
    emptydrops = [reverse_complement(i) for i in emptydrops.values]
    return set(emptydrops)

def bam_to_fastq(read):
    name = read.query_name
    seq = read.query_sequence
    qual = read.qual
    return "@{}\n{}\n+\n{}\n".format(name,seq,qual)

def read_fltnc(sample,emptydrops):
    
    dic_bam_per_cell=defaultdict(lambda:[])
    dic_UMI_per_cell=defaultdict(lambda:defaultdict(lambda:True))

    bams = glob.glob('input_flntc/{}*.bam'.format(sample))
    
    for bamfile in bams:
        samfile = pysam.AlignmentFile(bamfile, "rb", check_sq=False, threads=args.cpu)
        for READ in samfile:
            try:
                BC = READ.get_tag('XC')
                UMI = READ.get_tag('XM')
            except KeyError:
                continue
            if str(BC) in emptydrops:
                READ.query_name = READ.query_name + '_' + BC
                if dic_UMI_per_cell[BC][UMI]:
                    dic_bam_per_cell[BC].append(READ)
                    dic_UMI_per_cell[BC][UMI]=False

    print("{} mean reads per empty droplet: {} reads".format(sample,
                                                            np.mean([len(i) for i in dic_bam_per_cell.values()])))
    print("{} # dead cell: {}".format(sample, len(dic_bam_per_cell.values())))
    print()
    
    return dic_bam_per_cell
    
                
def main(args):

    sample = args.sample

    emptydrops = get_emptydrops(args.emptydroplets)

    dic_bam_per_cell= read_fltnc(sample,emptydrops)   

    path = 'emptydroplets/{}.fastq.gz'.format(sample)
    with gzip.open(path, 'wt') as fastq:
        for BC in dic_bam_per_cell:
            for read in dic_bam_per_cell[BC]:
                fastq.write(bam_to_fastq(read))


def parse_args():
    parser = argparse.ArgumentParser(
        prog='empty_droplets_listing.py', 
        usage='python3 empty_droplets_listing.py --emptydroplets <emptydrops.tsv> --sample <sample> ',
        description='Divides bamfiles per cell prior to UMI deduplication'
    )
    parser.add_argument(
        '--emptydroplets', type=str,
        help='path(s) to directory containing barcodes of cellranger empty droplets, tsv format'
    )
    parser.add_argument(
        '--sample', type=str,
        help='sample name (should not contain ".")'
    )
    parser.add_argument(
        '--cpu', type=int, default=8,
        help='# CPUs to use'
    )


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
