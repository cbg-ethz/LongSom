import pandas as pd
import numpy as np
import pysam
from pathlib import Path
import glob
import argparse
import gzip
from collections import defaultdict
import re


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

def get_ReadBase(READ,mut_pos):

    SEQ = READ.query_sequence
    START = READ.reference_start
    CIGAR = READ.cigarstring
    if CIGAR is None:
        return 'NA','NA'
    NAME = READ.query_name

    CELL = NAME.split('_')[-1]

    pos_genome = START #initialize position on genome
    if pos_genome+1>mut_pos:
        return 'NA','NA'
    pos_read = 0 #initialize position on read
    
    cigartuples = get_cigartuple(CIGAR)
    
    pos_read = get_readpos(cigartuples,pos_read,pos_genome,mut_pos)
    
    if pos_read == -1:
        return 'NA','NA'
    
    base = SEQ[pos_read-1]
    return CELL, base

def count_bases(sam,idx,mchr,mut_pos,ref_base,alt_base,cellinfo):
    
    for READ in sam.fetch(mchr,mut_pos,mut_pos+1):
        CELL, base = get_ReadBase(READ,mut_pos)
        if base =='NA':
            continue
        elif base == ref_base:
            cellinfo[idx][CELL]['REF']+=1
        elif base == alt_base:
            cellinfo[idx][CELL]['ALT']+=1
        else:
            cellinfo[idx][CELL]['OTHER']+=1
        cellinfo[idx][CELL]['TOT']+=1

    return cellinfo

def muts_per_cells(vcf, sam):
    cellinfo = {}

    for idx,line in vcf.iterrows():
        mchr = line['#CHROM']
        mut_pos = line['POS']
        ref_base = line['REF']
        alt_base = line['ALT']
        idx = line['INDEX']
        cellinfo[idx] = defaultdict(lambda:{'REF':0,'ALT':0,'OTHER':0,'TOT':0})

        cellinfo = count_bases(sam,idx,mchr,mut_pos,ref_base,alt_base,cellinfo)

    return cellinfo

def main(args):

    sam = pysam.AlignmentFile(args.bam, "rb", threads=args.cpu)

    stats = sam.get_index_statistics()
    print(args.sample)
    print("Reads mapping to chrM: {}".format(stats[324].mapped))
    print("Total reads mapped: {}".format(sam.mapped))
    print("Perc. reads ampping to chrM: {}".format(stats[324].mapped/sam.mapped))

    vcf = pd.read_csv(args.vcf, sep = '\t')

    cellinfo = muts_per_cells(vcf, sam)
    d = {'Barcodes':[],
        'CellType':[],
        'ALT':[],
        'VAF':[],
        'INDEX':[],}

    for IDX in cellinfo:
        print("Covered cells for {} : {}".format(IDX, len(cellinfo[IDX])))
        print("Mean ALT frac for {} : {}".format(IDX, 
            np.sum([cellinfo[IDX][c]['ALT'] for c in cellinfo[IDX]]) / np.sum([cellinfo[IDX][c]['TOT'] for c in cellinfo[IDX]])
                                                 ))
        print("Mean REF frac for {} : {}".format(IDX, 
            np.sum([cellinfo[IDX][c]['REF'] for c in cellinfo[IDX]]) / np.sum([cellinfo[IDX][c]['TOT'] for c in cellinfo[IDX]])
                                                 ))
        
        cells = cellinfo[IDX].keys()
        alts,tots,fracs = [],[],[]

        for CELL in cellinfo[IDX]:
            alt = cellinfo[IDX][CELL]['ALT']
            alts.append(alt)
            tot = cellinfo[IDX][CELL]['TOT']
            tots.append(tot)
            frac = alt / tot
            fracs.append(frac)

        if 'Om' in args.sample:
            ctype = 'Empty droplets\n Normal'
        elif 'Tum' in args.sample:
            ctype = 'Empty droplets\n Tumor'



        df = pd.DataFrame({'Barcodes':cells, 'ALT'.format(IDX):alts, 'VAF':fracs})
        df['VAF'].round(2)
        df['CellType'] = ctype
        df['INDEX'] = IDX
        df.to_csv("results/{}/{}.{}.emptydrops.csv".format(args.sample, args.sample, IDX), index = False)

        d['Barcodes'] += list(cells)
        d['CellType'] += [ctype]*len(list(cells))
        d['ALT'] += list(alts)
        d['VAF'] += list(fracs)
        d['INDEX'] += [IDX]*len(list(cells))

    df = pd.DataFrame(d)
    df['VAF'].round(2)
    df.to_csv("results/{}/{}.emptydrops.tsv".format(args.sample, args.sample), sep ='\t', index = False)


def parse_args():
    parser = argparse.ArgumentParser(
        prog='split_cells_bam.py', 
        usage='python3 split_cells_bam.py --bc_dir <bc_dir> --sample <sample> ',
        description='Divides bamfiles per cell prior to UMI deduplication'
    )
    parser.add_argument(
        '--bam', type=str,
        help='indexed and sorted input bam'
    )
    parser.add_argument(
        '--vcf', type=str,
        help='LongSom output muts vcf'
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
    



        