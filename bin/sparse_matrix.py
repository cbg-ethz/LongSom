from collections import defaultdict, Counter
import pandas as pd
import numpy as np
import argparse
import pysam
import re
from itertools import compress



def get_cigartuple(CIGAR):
    l=re.split('(\d+)', CIGAR)
    cigartuples=[[int(l[i+1]),l[i+2]] for i in range(0,len(l)-2,2)]
    return cigartuples


def get_readpos(cigartuples,pos_read,pos_genome,mut_pos):
    
    for i in cigartuples:
        if i[1]=='H':
            continue
        elif i[1]=='S':
            pos_read += i[0]
        elif i[1]=='M':
            pos_genome += i[0] #progress on genome
            if pos_genome >= mut_pos: #if exon exceed mutation pos, find mutation pos on read
                pos_read += i[0] - (pos_genome - mut_pos) #mutation pos on read = actual pos + distance to pos in exon 
                return pos_read
            pos_read += i[0]    
        elif i[1]=='N':
            pos_genome += i[0]
            if pos_genome >= mut_pos:
                return -1
        elif i[1]=='I':
            pos_read += i[0]
        elif i[1]=='D':
            pos_genome += i[0]
            if pos_genome >= mut_pos:
                return -1

def get_celltype(CELL, celltypes):
    try:
        ctype = celltypes[CELL]
    except KeyError:
        ctype = 'NA'
    return ctype

def get_ReadBase_LR(READ,mut_pos,celltypes):

    SEQ = READ.query_sequence
    START = READ.reference_start
    CIGAR = READ.cigarstring
    if CIGAR is None:
        return 'NA','NA','NA'

    CELL = READ.get_tag('CB')
    
    CTYPE = get_celltype(CELL, celltypes)
    if CTYPE == 'NA':
        return 'NA','NA','NA'

    pos_genome = START #initialize position on genome

    # fetch(mchr,mut_pos,mut_pos+1) also captures reads starting 1bp after SNV locus
    # fetch(mchr,mut_pos,mut_pos) does not work
    # this is my current workaround to filter reads starting after locus:
    if pos_genome+1>mut_pos: 
        return 'NA','NA','NA'

    pos_read = 0 #initialize position on read
    cigartuples = get_cigartuple(CIGAR)
    pos_read = get_readpos(cigartuples,pos_read,pos_genome,mut_pos)
    
    if pos_read == -1:
        return 'NA','NA','NA'
    
    base = SEQ[pos_read-1]
    return CTYPE, CELL, base


def count_bases(sam_tum,mchr,mut_pos,celltypes,alt_base,mutcells, covcells,idx):
    
    for READ in sam_tum.fetch(mchr,mut_pos,mut_pos+1):
        CTYPE, CELL, base = get_ReadBase_LR(READ,mut_pos,celltypes)

        #ignore reads from cells not in list
        if CTYPE == 'NA':
            continue

        covcells[CELL][idx]+=1
        if base == alt_base:
            mutcells[CELL][idx]+=1

    return mutcells, covcells

def div_cov(mutcells, covcells):
    cov = pd.DataFrame(covcells)
    mut = pd.DataFrame(mutcells)
    freq = (mut/cov).replace(np.inf, -1)
    freq = freq.replace(np.nan, -1)
    return freq, cov, mut


def main(args):

    names = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SPL']
    muts = pd.read_csv(args.tum_vcf, comment='#', delim_whitespace=True,
             header=None, names=names, compression='gzip', encoding = "ISO-8859-1")

    ### REMOVING INDELS, HIGHLY RECOMMENDED FOR LR
    if args.snv:
        mask = (muts ['ALT'].str.len() == 1) & (muts ['REF'].str.len() == 1)
        muts  = muts .loc[mask]

    ### LOADING CELLTYPES
    bc_to_pheno = pd.read_csv(args.ctypes, sep='\t')
    barcodes = list(bc_to_pheno['barcodes'])
    celltypes = dict(zip(bc_to_pheno['barcodes'],bc_to_pheno['celltype_final']))


    ### PREPARING EMPTY SPARSE MATRICES (one coverage one mutated)
    muts['INDEX'] = muts['#CHROM'] + ':' + muts['POS'].astype(str) + ':' + muts['ALT'] 
    muts.set_index('INDEX', inplace=True)
    covcells_LR = {j:{i:0 for i in muts.index} for j in barcodes}
    mutcells_LR = {j:{i:0 for i in muts.index} for j in barcodes}

    ###LOADING BAM
    sam_tum_LR = pysam.AlignmentFile(args.bam_tum_LR, "rb", threads=args.cpu)

    ### Iterating through all SNV entries from VCF
    for index,mut in muts.iterrows():

        alt_base = mut['ALT']
        mchr = mut['#CHROM'] 
        mut_pos = mut['POS']

        #Tum LR scRNA
        mutcells_LR, covcells_LR =count_bases(
            sam_tum_LR,mchr,mut_pos,celltypes,alt_base, 
            mutcells_LR, covcells_LR,index)

    freqtum, covtum, muttum = div_cov(mutcells_LR, covcells_LR)

    freqtum.to_csv('VAF_sparse_matrix.tsv',sep='\t')
    covtum.to_csv('Coverage_sparse_matrix.tsv',sep='\t')
    muttum.to_csv('Variant_sparse_matrix.tsv',sep='\t')
    
def parse_args():
    parser = argparse.ArgumentParser(
        prog='filter_dist_nontum_final.py', 
        usage='python3 FMI_mutations_to_cells.py --bam <sorted.bam> --csv FMI_mutations.csv',
        description='Creates csv files of mutated cell names and their mutations'
    )
    parser.add_argument(
        '--tum_vcf', type=str, help='tumor cells vcf.gz file'
    )
    parser.add_argument(
        '--snv', action='store_true', help='snv only'
    )
    parser.add_argument(
        '--bam_tum_LR', type=str,
        help='Absolute or relative path(s) to tumor bam file'
    )
    parser.add_argument(
        '--ctypes', type=str,
        help='Absolute or relative path(s) to input barcode-to-ctype file'
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
