from collections import defaultdict, Counter
import pandas as pd
import numpy as np
import argparse
import pysam
import re
import time
from itertools import compress
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from gnomad_db.database import gnomAD_DB

def none_or_str(value):
    if value == 'None':
        return None
    else:
        return value

def gnomad(line):
    if line['#CHROM'] == 'chrM':
        return 0
    try: 
        l = float(line['INFO'].split("gnomad_AF=")[1].split(';')[0])
    except IndexError:
        l = 0
    return l

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq)
    letters = [complement[base] for base in bases] 
    letters = ''.join(letters)
    reverse = letters[::-1]
    return reverse

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
    NAME = READ.query_name

    CELL = NAME.split('.')[-1]
    CELL = reverse_complement(CELL)
    
    CTYPE = get_celltype(CELL, celltypes)
    if CTYPE == 'NA':
        return 'NA','NA','NA'

    pos_genome = START #initialize position on genome
    if pos_genome+1>mut_pos:
        return 'NA','NA','NA'
    pos_read = 0 #initialize position on read
    
    cigartuples = get_cigartuple(CIGAR)
    
    pos_read = get_readpos(cigartuples,pos_read,pos_genome,mut_pos)
    
    if pos_read == -1:
        return 'NA','NA','NA'
    
    base = SEQ[pos_read-1]
    return CTYPE, CELL, base

def get_ReadBase_SR(READ,mut_pos,celltypes):

    SEQ = READ.query_sequence
    START = READ.reference_start
    CIGAR = READ.cigarstring
    if CIGAR is None:
        return 'NA','NA','NA'
    NAME = READ.query_name

    try:
        CELL = READ.get_tag('CB').split('-')[0]
    except KeyError:
        return 'NA','NA','NA'

    CTYPE = get_celltype(CELL, celltypes)
    if CTYPE == 'NA':
        return 'NA','NA','NA'

    pos_genome = START #initialize position on genome
    if pos_genome+1>mut_pos:
        return 'NA','NA','NA'
    pos_read = 0 #initialize position on read
    
    cigartuples = get_cigartuple(CIGAR)
    
    pos_read = get_readpos(cigartuples,pos_read,pos_genome,mut_pos)
    
    if pos_read == -1:
        return 'NA','NA','NA'
    
    base = SEQ[pos_read-1]
    return CTYPE, CELL, base

def get_ReadBase_scDNA(READ,mut_pos,celltypes):

    SEQ = READ.query_sequence
    START = READ.reference_start
    CIGAR = READ.cigarstring
    if CIGAR is None:
        return 'NA','NA','NA','NA'
    NAME = READ.query_name

    try:
        CELL = READ.get_tag('CB')
    except KeyError:
        return 'NA','NA','NA','NA'

    CTYPE = get_celltype(CELL, celltypes)
    if CTYPE == 'NA':
        return 'NA','NA','NA', 'NA'

    pos_genome = START #initialize position on genome
    if pos_genome+1>mut_pos:
        return 'NA','NA','NA', 'NA'
    pos_read = 0 #initialize position on read
    
    cigartuples = get_cigartuple(CIGAR)
    
    pos_read = get_readpos(cigartuples,pos_read,pos_genome,mut_pos)
    
    if pos_read == -1:
        return 'NA','NA','NA','NA'
    
    base = SEQ[pos_read-1]
    return NAME, CTYPE, CELL, base

def count_bases(sam_tum,mchr,mut_pos,celltypes,ref_base,alt_base,cell_info,mutcells, covcells, biopsy, readtype):
    
    for READ in sam_tum.fetch(mchr,mut_pos,mut_pos+1):
        if readtype == 'LR':
            CTYPE, CELL, base = get_ReadBase_LR(READ,mut_pos,celltypes)
            idx = ':'.join([mchr,str(mut_pos),alt_base])
        elif readtype == 'SR':
            if mchr == 'MT':
                mchr = 'M'
            CTYPE, CELL, base = get_ReadBase_SR(READ,mut_pos,celltypes)
            idx = ':'.join(['chr{}'.format(str(mchr)),str(mut_pos),alt_base])

        if CTYPE == 'NA':
            continue

        if biopsy == 'Tumor':
            if CTYPE=='HGSOC':
                covcells['Tumor'][CELL][idx]+=1
                cell_info['Tumor_{}'.format(readtype)]['DP']+=1
                if base == alt_base:
                    mutcells['Tumor'][CELL][idx]+=1
                    cell_info['Tumor_{}'.format(readtype)]['ALT']+=1
                elif base == ref_base:
                    cell_info['Tumor_{}'.format(readtype)]['REF']+=1

            else:
                covcells['NonTumor'][CELL][idx]+=1
                cell_info['NonTumor_{}'.format(readtype)]['DP']+=1
                if base == alt_base:
                    cell_info['NonTumor_{}'.format(readtype)]['ALT']+=1
                    mutcells['NonTumor'][CELL][idx]+=1
                elif base == ref_base:
                    cell_info['NonTumor_{}'.format(readtype)]['REF']+=1

        else:
            covcells['Distal'][CELL][idx]+=1
            cell_info['Distal_{}'.format(readtype)]['DP']+=1
            if base == alt_base:
                cell_info['Distal_{}'.format(readtype)]['ALT']+=1
                mutcells['Distal'][CELL][idx]+=1
            elif base == ref_base:
                cell_info['Distal_{}'.format(readtype)]['REF']+=1

    return cell_info, mutcells, covcells



def count_bases_scDNA(sam_tum,mchr,mut_pos,celltypes,ref_base,alt_base,cell_info):
    reads = []
    for READ in sam_tum.fetch(mchr,mut_pos,mut_pos+1):
        NAME, CTYPE, CELL, base = get_ReadBase_scDNA(READ,mut_pos,celltypes)

        if NAME in reads:
            continue
        reads.append(NAME)

        if CTYPE == 'NA':
            continue

        cell_info[CTYPE]['DP']+=1
        if base == alt_base:
            cell_info[CTYPE]['ALT']+=1
        elif base == ref_base:
            cell_info[CTYPE]['REF']+=1

    return cell_info


def div_cov(mutcells, covcells, smpl):
    cov = pd.DataFrame(covcells[smpl])
    mut = pd.DataFrame(mutcells[smpl])
    ncells = len(cov.columns) 
    freq = (mut/cov).replace(np.inf, 'NoCov')
    freq = freq.replace(np.nan, 'NoCov')
    freq['MutatedCells_{}'.format(smpl)] = freq.replace({'NoCov':0}).gt(.3).sum(axis=1)
    freq['FracMut_{}'.format(smpl)] = freq['MutatedCells_{}'.format(smpl)]/ncells
    freq['FracCov_{}'.format(smpl)] = (ncells - (freq == 'NoCov').sum(axis=1))/ncells
    freq['FracMutCov_{}'.format(smpl)] = freq['FracMut_{}'.format(smpl)]/freq['FracCov_{}'.format(smpl)]
    freq = freq.replace(np.inf,0)
    freq = freq.replace(np.nan,0)
    return freq, cov, mut

def add_fusions(freq,fusions_file):
    freq = freq.replace({'NoCov':np.nan})

    fusions =  pd.read_csv(fusions_file, sep = '\t')

    d = defaultdict(lambda: [])
    for idx, row in fusions.iterrows():
        d[row['FusionName']].append(row['barcodes'])
        
    nd = {}
    for k in d:
        nd[k] = list(set(d[k]))

    f = {j:{i:0 for i in nd} for j in freq.columns}

    for k in nd:
        for v in nd[k]:
            if v in freq.columns:
                f[v][k]+=1
    fusion = pd.DataFrame(f)
    fusion = fusion.replace({0:np.nan})
    fusion['Diff'] = -3
    freq = pd.concat([freq,fusion])

    return freq, fusion

def main(args):
    start = time.time()

    if args.blacklist is not None:
        blacklist = list(pd.read_csv(args.blacklist, sep = '\t', header=None, names=['BL'])['BL'])
    else:
        blacklist = []

    names = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SPL']
    tum_vcf = pd.read_csv(args.tum_vcf, comment='#', delim_whitespace=True,
             header=None, names=names, compression='gzip', encoding = "ISO-8859-1")

    if args.snv:
        mask = (tum_vcf['ALT'].str.len() == 1) & (tum_vcf['REF'].str.len() == 1)
        tum_vcf = tum_vcf.loc[mask]
        
    keys = ['#CHROM', 'POS', 'ALT']
    i2 = tum_vcf.set_index(keys).index

    if args.dist_vcf is not None:
        dist_vcf = pd.read_csv(args.dist_vcf, comment='#', delim_whitespace=True,
         header=None, names=names, compression='gzip', encoding = "ISO-8859-1")
        i1 = dist_vcf.set_index(keys).index
        tum_vcf = tum_vcf[~i2.isin(i1)]

    tum_vcf = tum_vcf[~tum_vcf['INFO'].str.contains("RNAEDIT")]
    tum_vcf[['GT','AD','DP','GQ','PL']] = tum_vcf['SPL'].str.split(':',expand=True)
    tum_vcf[['refDP','altDP']] = tum_vcf['AD'].str.split(',',expand=True)
    tum_vcf['DP'] = tum_vcf['DP'].astype('int')
    tum_vcf = tum_vcf[tum_vcf['DP']>args.mindepth]
    tum_vcf['altDP'] = tum_vcf['altDP'].astype('int')
    tum_vcf['refDP'] = tum_vcf['refDP'].astype('int')
    tum_vcf = tum_vcf[tum_vcf['refDP']>=args.minref]
    tum_vcf = tum_vcf[tum_vcf['altDP']>=args.minalt]
    tum_vcf = tum_vcf[names]

    muts = tum_vcf
    muts['INDEX'] = muts['#CHROM'] + ':' + muts['POS'].astype(str) + ':' + muts['ALT']

    ###germline filtering with gnomAD
    muts = muts.reset_index(drop=True)
    database_location = args.gnomAD
    db = gnomAD_DB(database_location, gnomad_version="v4")
    gnom = muts[['#CHROM','POS','REF','ALT']]
    gnom.columns = ['chrom','pos','ref','alt']
    muts['gnomAD'] = db.get_info_from_df(gnom, "AF")
    muts['gnomAD'] = muts['gnomAD'].replace(np.nan,0)
    muts = muts[muts['gnomAD']<0.01]
    ###
    
    ### removing mutations close to eachother
    idx = muts['INDEX']
    a=[]
    for i in idx:
        chr,pos,base=i.split(':')
        a.append((chr,pos,base))

    b = sorted(a, key=lambda x: (x[0],x[1]))

    trash=[]

    for (chr,pos,base), (chr2,pos2,base2) in zip(b, b[1:]):
        if chr==chr2:
            if chr == 'chrM':
                continue
            if abs(int(pos)-int(pos2))<10000:
                trash.append(':'.join([chr,pos,base]))
                trash.append(':'.join([chr2,pos2,base2]))

    muts = muts[~muts['INDEX'].isin(trash)]
    ###
    
    ###LOADING BAMS
    sam_tum_LR = pysam.AlignmentFile(args.bam_tum_LR, "rb", threads=args.cpu)
    if args.bam_dist_LR is not None:
        sam_dist_LR = pysam.AlignmentFile(args.bam_dist_LR, "rb", threads=args.cpu)
    sam_tum_SR = pysam.AlignmentFile(args.bam_tum_SR, "rb", threads=args.cpu)
    if args.bam_dist_SR is not None:
        sam_dist_SR = pysam.AlignmentFile(args.bam_dist_SR, "rb", threads=args.cpu)
    sam_scDNA = pysam.AlignmentFile(args.bam_scDNA, "rb", threads=args.cpu)
    ###

    ctypes_categ = ["HGSOC", "Mesothelial.cells", "Fibroblasts", "T.NK.cells", "Myeloid.cells", "B.cells", "Endothelial.cells"]

    bc_to_pheno = pd.read_csv(args.ctypes, sep='\t')
    bc_to_pheno['celltype_final'] = pd.Categorical(bc_to_pheno['celltype_final'], ctypes_categ)
    bc_to_pheno = bc_to_pheno.sort_values(by = ['biopsy', 'celltype_final'])
    barcodes = list(bc_to_pheno['barcodes'])

    cells_tum = bc_to_pheno[bc_to_pheno['celltype_final']=='HGSOC']
    cells_nontum = bc_to_pheno[(bc_to_pheno['celltype_final']!='HGSOC') & (bc_to_pheno['biopsy']=='Tumor')]
    cells_dist = bc_to_pheno[bc_to_pheno['biopsy']=='Distal']
    
    cells = {'Tumor': list(cells_tum['barcodes']),
            'NonTumor':list(cells_nontum['barcodes']),
            'Distal': list(cells_dist['barcodes'])}

    celltypes = {'Tumor': dict(zip(cells_tum['barcodes'],cells_tum['celltype_final'])),
            'NonTumor':dict(zip(cells_nontum['barcodes'],cells_nontum['celltype_final'])),
            'Distal': dict(zip(cells_dist['barcodes'],cells_dist['celltype_final']))}

    celltypes['Metastasis'] = celltypes['Tumor']|celltypes['NonTumor']


    covcells_LR = {'Tumor':{j:{i:0 for i in muts['INDEX']} for j in cells['Tumor']},
                'NonTumor':{j:{i:0 for i in muts['INDEX']} for j in cells['NonTumor']},
                'Distal':{j:{i:0 for i in muts['INDEX']} for j in cells['Distal']}}

    mutcells_LR = {'Tumor':{j:{i:0 for i in muts['INDEX']} for j in cells['Tumor']},
                'NonTumor':{j:{i:0 for i in muts['INDEX']} for j in cells['NonTumor']},
                'Distal':{j:{i:0 for i in muts['INDEX']} for j in cells['Distal']}}

    covcells_SR = {'Tumor':{j:{i:0 for i in muts['INDEX']} for j in cells['Tumor']},
                'NonTumor':{j:{i:0 for i in muts['INDEX']} for j in cells['NonTumor']},
                'Distal':{j:{i:0 for i in muts['INDEX']} for j in cells['Distal']}}

    mutcells_SR = {'Tumor':{j:{i:0 for i in muts['INDEX']} for j in cells['Tumor']},
                'NonTumor':{j:{i:0 for i in muts['INDEX']} for j in cells['NonTumor']},
                'Distal':{j:{i:0 for i in muts['INDEX']} for j in cells['Distal']}}

    clones = pd.read_csv(args.scDNA_clones, sep='\t')
    clones = dict(zip(clones['Barcode'],clones['Clone']))
    clone_list = list(set(clones.values()))

    t1 = time.time()
    print('Time to load bams: ',t1-start)

    for index,mut in muts.iterrows():

        cell_info = {'Tumor_LR':{'DP':0,'REF':0, 'ALT':0},
                        'NonTumor_LR':{'DP':0,'REF':0, 'ALT':0},
                        'Distal_LR':{'DP':0,'REF':0, 'ALT':0},
                        'Tumor_SR':{'DP':0,'REF':0, 'ALT':0},
                        'NonTumor_SR':{'DP':0,'REF':0, 'ALT':0},
                        'Distal_SR':{'DP':0,'REF':0, 'ALT':0}
                        }

        for clone in clone_list:
            cell_info[clone] = {'DP':0,'REF':0, 'ALT':0}

        cell_info['scDNA_Tumor'] = {'DP':0,'REF':0, 'ALT':0}
        cell_info['scDNA_NonTumor'] = {'DP':0,'REF':0, 'ALT':0}

        alt_base = mut['ALT']
        ref_base = mut['REF']
        mchr = mut['#CHROM'] 
        mut_pos = mut['POS']

        #Tum LR scRNA
        cell_info, mutcells_LR, covcells_LR =count_bases(
            sam_tum_LR,mchr,mut_pos,celltypes['Metastasis'],ref_base,alt_base,cell_info, 
            mutcells_LR, covcells_LR,'Tumor', 'LR')

        #Distal LR scRNA
        if args.bam_dist_LR is not None:
            cell_info, mutcells_LR, covcells_LR = count_bases(
                sam_dist_LR,mchr,mut_pos,celltypes['Distal'],ref_base,alt_base,cell_info,
                 mutcells_LR, covcells_LR ,'Distal', 'LR')

        #Tum SR scRNA
        mchr_SR = re.split('chr', mchr)[1]
        if mchr_SR == 'M':
            mchr_SR = 'MT'
        cell_info, mutcells_SR, covcells_SR =count_bases(
            sam_tum_SR,mchr_SR,mut_pos,celltypes['Metastasis'],ref_base,alt_base,cell_info,
            mutcells_SR, covcells_SR, 'Tumor','SR')
        #Distal SR scRNA
        if args.bam_dist_SR is not None:
            cell_info, mutcells_SR, covcells_SR = count_bases(
                sam_dist_SR,mchr_SR,mut_pos,celltypes['Distal'],ref_base,alt_base,cell_info,
                 mutcells_SR, covcells_SR,'Distal', 'SR')

           #Tum SR scDNA
        cell_info =count_bases_scDNA(
            sam_scDNA,mchr,mut_pos,clones,ref_base,alt_base,cell_info)
        
        
        for clone in clone_list:
            if '_Tum' in clone:
                cell_info['scDNA_Tumor']['DP']+=cell_info[clone]['DP']
                cell_info['scDNA_Tumor']['REF']+=cell_info[clone]['REF']
                cell_info['scDNA_Tumor']['ALT']+=cell_info[clone]['ALT']
            elif '_NonTum' in clone:
                cell_info['scDNA_NonTumor']['DP']+=cell_info[clone]['DP']
                cell_info['scDNA_NonTumor']['REF']+=cell_info[clone]['REF']
                cell_info['scDNA_NonTumor']['ALT']+=cell_info[clone]['ALT']

        for cond in cell_info:
            muts.loc[index,cond] = '{}:{},{}'.format(cell_info[cond]['DP'], cell_info[cond]['REF'], cell_info[cond]['ALT'])

        if cell_info['scDNA_NonTumor']['DP'] == 0:
            muts.loc[index,'scDNASupport_NonTumor'] = np.nan
        else:
            muts.loc[index,'scDNASupport_NonTumor'] = cell_info['scDNA_NonTumor']['ALT']/cell_info['scDNA_NonTumor']['DP']
        
        if cell_info['scDNA_Tumor']['DP'] == 0:
            muts.loc[index,'scDNASupport_Tumor'] = np.nan
        else:
            muts.loc[index,'scDNASupport_Tumor'] = cell_info['scDNA_Tumor']['ALT']/cell_info['scDNA_Tumor']['DP']

    
    t2=time.time()
    print('Time to read bams:',t2-t1)       
    muts.to_csv('{}.barcodes_added.vcf'.format(args.sample), sep = '\t', index = False)
    print('#Mutations after filtering SNVs only and min alt/ref: ', len(muts))

    freqtum, covtum, muttum = div_cov(mutcells_LR, covcells_LR, 'Tumor')
    freqnontum, covnontum, mutnontum = div_cov(mutcells_LR, covcells_LR, 'NonTumor')
    freqdist, covdist, mutdist = div_cov(mutcells_LR, covcells_LR, 'Distal')

    freqtum_SR,_,_ = div_cov(mutcells_SR, covcells_SR, 'Tumor')
    freqnontum_SR,_,_ = div_cov(mutcells_SR, covcells_SR, 'NonTumor')
    freqdist_SR,_,_ = div_cov(mutcells_SR, covcells_SR, 'Distal')

    freq_SR = pd.concat([freqtum_SR, freqnontum_SR, freqdist_SR], axis = 1)

    freqtum.to_csv('{}.tum.{}muts.tsv'.format(args.sample, str(args.maxmuts)), sep = '\t', index=True)
    freqnontum.to_csv('{}.nontum.{}muts.tsv'.format(args.sample, str(args.maxmuts)), sep = '\t', index=True)
    freqdist.to_csv('{}.dist.{}muts.tsv'.format(args.sample, str(args.maxmuts)), sep = '\t', index=True)
    
    freq = pd.concat([freqtum, freqnontum, freqdist], axis = 1)

    cols = [i for i in list(freq.columns) if i in barcodes]
    n = len(cols)

    freq, fusion = add_fusions(freq,args.fusions)

    ### HIGH CONFIDENCE CANCER VARIANTS IDENTIFICATION

    if args.bam_dist_LR is not None:
        tum_muts = list(freq[(freq['FracMut_Tumor']>0.05) & 
                (freq['FracMutCov_Tumor']>0.2) &
                (freq['FracMutCov_NonTumor']<0.05) &
                (freq['FracCov_NonTumor']>0.01) &
                (freq['FracMutCov_Distal']<0.001)].index)+list(fusion.index)
    else:
        tum_muts = list(freq[(freq['FracMut_Tumor']>0.05) & 
                (freq['FracMutCov_Tumor']>0.2) &
                (freq['FracMutCov_NonTumor']<0.05) &
                (freq['FracCov_NonTumor']>0.01)].index)+list(fusion.index)

    freq = freq[cols]
    freq = freq.replace({'NoCov':np.nan})

    #redefinition of cell typing by mutation per cell typed as tumor
    df = freq.loc[tum_muts]
    df.replace(3, np.nan, inplace=True)
    df[df <0.3] = 0
    df[df >=0.3] = 1
    df.loc['Sum',:] = df.sum(axis=0)
    df.loc['Covered',:] = df.count()-1
    df.loc['Perc'] = df.loc['Sum'] / df.loc['Covered']
    tumor_cells = list(df.loc[:,(df.loc['Sum']>=2)].columns)
    ###retrieving tumor cells using more permissive filter:
    '''
    for i in list(cells_tum['barcodes']):
        if df.loc['Sum',i]>=1:
            if i not in tumor_cells:
                tumor_cells.append(i)
    '''

    #tumor_cells = df.loc[:,(df.loc['Perc']>.2) & (df.loc['Covered']>(len(df.index)-3)*0.05)].columns
    #tumor_cells = df.loc[:,(df.loc['Perc']>.5) & (df.loc['Sum']>(len(df.index)-3)*0.05)].columns
    #50% of covered loci are mutated and 5% of all tumor loci covered (previously mutated)
    for i in tumor_cells:
        print(i)

    nontumor_cells = [i for i in cols if i not in tumor_cells]

    #calculation of tumor and non tumor mutation fraction based on new celltyping
    freq['Mutcells_Tumor'] = freq[tumor_cells].gt(.3,axis = 0).sum(axis = 1)
    freq['FracMut_Tumor'] = freq[tumor_cells].gt(.3,axis = 0).sum(axis = 1)/len(tumor_cells)
    freq['FracCov_Tumor'] = freq[tumor_cells].count(axis=1)/len(tumor_cells)
    freq['FracMutCov_Tumor'] = freq['FracMut_Tumor']/freq['FracCov_Tumor']
    freq['FracMut_NonTumor'] = freq[nontumor_cells].gt(.3,axis = 0).sum(axis=1)/len(nontumor_cells)
    freq['FracCov_NonTumor'] = (freq[nontumor_cells].count(axis=1))/len(nontumor_cells)
    freq['FracMutCov_NonTumor'] = freq['FracMut_NonTumor']/freq['FracCov_NonTumor']
    freq = freq.replace(np.inf,0)
    freq_chrM = freq.loc[[i for i in freq.index if 'chrM' in i]]
    freq_chrM['Diff'] = abs(freq_chrM['FracMutCov_NonTumor']-freq_chrM['FracMutCov_Tumor'])
    keep_chrm =list(freq_chrM[freq_chrM['Diff']>0.2].index)
    
    # #fus_save = freq.loc[fusion.index]
    # tum_muts = list(freq[(freq['FracMut_Tumor']>0.05) & 
    #             (freq['FracMutCov_NonTumor']<0.01) &
    #             (freq['FracCov_NonTumor']>0.005)].index)+list(fusion.index)


    # save=list(set(tum_muts+keep_chrm))

    # #pre-filter to compute real frac_tumor
    # freq = freq.loc[save]
    # #remove empty cells
    # freq = freq.loc[:, freq.count()>=5]
    # #freq = freq.loc[:, (freq != 3).sum(axis=0)>=3]


    # #recompute fractions after removing empty cells
    # tumor_cells = [i for i in tumor_cells if i in freq.columns]
    # nontumor_cells = [i for i in nontumor_cells if i in freq.columns]

    # freq['FracMut_Tumor'] = freq[tumor_cells].gt(.3,axis = 0).sum(axis = 1)/len(tumor_cells)
    # freq['FracCov_Tumor'] = freq[tumor_cells].count(axis=1)/len(tumor_cells)
    # freq['FracMutCov_Tumor'] = freq['FracMut_Tumor']/freq['FracCov_Tumor']
    # freq['FracMut_NonTumor'] = freq[nontumor_cells].gt(.3,axis = 0).sum(axis=1)/len(nontumor_cells)
    # freq['FracCov_NonTumor'] = (freq[nontumor_cells].count(axis=1))/len(nontumor_cells)
    # freq['FracMutCov_NonTumor'] = freq['FracMut_NonTumor']/freq['FracCov_NonTumor']
    # freq = freq.replace(np.inf,0)
    
    #fus_save = freq.loc[fusion.index]
    tum_muts = list(freq[((freq['FracMut_Tumor']>0.05) |  (freq['Mutcells_Tumor']>=args.mincells)) & 
                (freq['FracMutCov_NonTumor']<0.01) &
                (freq['FracCov_NonTumor']>0.01)].index)+list(fusion.index)
    freq.to_csv('{}.freq.tsv'.format(args.sample), sep = '\t', index=True)

    save=list(set(tum_muts+keep_chrm))
    freq = freq.loc[save]



    freq['Diff'] =  abs(freq['FracMutCov_NonTumor'] - freq['FracMutCov_Tumor'])
    for i in fusion.index:
        freq.loc[i,'Diff'] = -3

    #mean mutfrac per ctype
    #freq, fusion = add_fusions(freq,args.fusions)
    #freq = pd.concat([freq,fusion])
    freq['VAF_Tumor'] = freq[tumor_cells].replace('NoCov',np.nan).apply(lambda row: np.nanmean(row), axis=1).replace(np.nan,0)
    freq['VAF_NonTumor'] = freq[nontumor_cells].replace('NoCov',np.nan).apply(lambda row: np.nanmean(row), axis=1).replace(np.nan,0)
    freq['VAF'] = freq.index +  '_' + ((freq['VAF_Tumor']*10).astype(int)/10).astype(str) + '_' + ((freq['VAF_NonTumor']*10).astype(int)/10).astype(str)
    VAF = dict(zip(freq.index,freq['VAF']))
    freq = freq.sort_values(by=['Diff'], ascending = [False])
    freq = freq[[i for i in freq.columns if i in cols]].head(args.maxmuts)
    freq = freq.loc[[locus for locus in freq.index if locus not in blacklist]] #filter blacklist loci
    freq = freq.loc[:, freq.count()>=3] #filter lowcov cells
    freq = freq.replace(np.nan,3)
    freq = freq[[i for i in freq if i not in cells['Distal']]]
    freq.rename(index=VAF).to_csv('{}.BnpC_input.{}muts.tsv'.format(args.sample, str(args.maxmuts)), sep = '\t',index=True)

    ###SR
    freq_SR = freq_SR[[i for i in freq.columns if i in freq_SR.columns]]
    freq_SR, _ = add_fusions(freq_SR, args.fusions_SR)
    freq_SR = freq_SR.loc[[i for i in freq.index if i in freq_SR.index]]
    freq_SR = freq_SR[~freq_SR.index.duplicated(keep='first')]
    for i in freq.index:
        if i not in freq_SR.index:
            freq_SR = freq_SR.reindex(freq_SR.index.values.tolist()+[i])
    freq_SR = freq_SR.reindex(freq.index)
    freq_SR = freq_SR.replace({np.nan:3})
    freq_SR.rename(index=VAF).to_csv('{}.BnpC_input_SR.{}muts.tsv'.format(args.sample, str(args.maxmuts)), sep = '\t',index=True)
    ###

    ###coverage and mutated reads only full matrices
    cov = pd.concat([covtum, covnontum, covdist], axis = 1)
    mut = pd.concat([muttum, mutnontum, mutdist], axis = 1)

    cov = cov.loc[[i for i in freq.index if i in cov.index]]
    mut = mut.loc[[i for i in freq.index if i in mut.index]]

    cov = cov[[i for i in freq.columns if i in cov.columns]]
    mut = mut[[i for i in freq.columns if i in mut.columns]]

    cov = cov.replace({np.nan:0})
    mut = mut.replace({np.nan:0})  

    cov.to_csv('{}.coverage.tsv'.format(args.sample), sep = '\t',index=True)
    mut.to_csv('{}.mutreadcount.tsv'.format(args.sample), sep = '\t',index=True)
    ###

    ### compute celltype correction data
    tumor_cells = [i for i in freq if i in tumor_cells]
    nontumor_cells = [i for i in freq if i in nontumor_cells]
    non_tumor_cells = [i for i in freq if ((i in nontumor_cells) & (i not in cells['Distal']))] #nontumor in Tumor biopsy only
    distal_cells = [i for i  in freq if i in cells['Distal']]

    bc_to_pheno_dict=dict(zip(bc_to_pheno['barcodes'],bc_to_pheno['celltype_final']))

    cells = tumor_cells + non_tumor_cells + distal_cells
    ctypes = ['Cancer']*len(tumor_cells) + ['Non-Cancer']*len(non_tumor_cells) + ['Normal']*len(distal_cells)


    bc_to_pheno_corr = dict(zip(cells,ctypes))
    bc_to_pheno['Corrected'] = bc_to_pheno['barcodes'].map(bc_to_pheno_corr)

    colors_corrected = {'Non-Cancer':'#94C773','Cancer':'#8F79A1','Normal':'grey'}

    bc_to_pheno['Colors_corrected'] = bc_to_pheno['Corrected'].map(colors_corrected)

    bc_to_pheno.to_csv('{}.celltype_correction.tsv'.format(args.sample), sep = '\t',index=False)

    print('cancer cells not HGSOC barcodes:')
    for i in tumor_cells:
        if bc_to_pheno_dict[i]!='HGSOC':
            print(i)
    print('')

    print('noncancer cells HGSOC barcodes:')
    for i in non_tumor_cells:
        if bc_to_pheno_dict[i]=='HGSOC':
            print(i)
    print('')

    tumor_celltypes=[bc_to_pheno_dict[i] for i in tumor_cells]
    nontumor_celltypes=[bc_to_pheno_dict[i] for i in non_tumor_cells]

    tumor_celltypes=Counter(tumor_celltypes)
    nontumor_celltypes=Counter(nontumor_celltypes)

    print('Perc TP',(tumor_celltypes['HGSOC']/tumor_celltypes.total())*100)
    print('Perc FN',((tumor_celltypes.total()-tumor_celltypes['HGSOC'])/tumor_celltypes.total())*100)

    print('Perc FP',(nontumor_celltypes['HGSOC']/nontumor_celltypes.total())*100)
    print('Perc TN',((nontumor_celltypes.total()-nontumor_celltypes['HGSOC'])/nontumor_celltypes.total())*100)

    ### Venn ploting
    print('entries covered LR: ',(freq != 3).sum().sum())
    print('entries covered SR: ',(freq_SR != 3).sum().sum())
    print('size: ', len(freq.columns)*len(freq.index))
    print('Perc entries covered LR: ',(freq != 3).sum().sum()/(len(freq.columns)*len(freq.index)))
    print('Perc entries covered SR: ',(freq_SR != 3).sum().sum()/(len(freq.columns)*len(freq.index)))

    freq[freq<3] = 1
    freq[freq==3] = 0

    freq_SR[freq_SR<3] = 2
    freq_SR[freq_SR==3] = 0

    venn = freq+freq_SR
    both = (venn==3).sum().sum()
    SR = (venn==2).sum().sum()
    LR = (venn==1).sum().sum()
    Nocov = (venn==0).sum().sum()

    print()
    print('Venn Diagram Values')
    print('Both :', both)
    print('LR only:', LR)
    print('SR only:', SR)
    print('Non covered:', Nocov)
    print('Perc cov LR : ', (LR+both)/(LR+both+SR+Nocov))
    print('Perc cov SR : ', (SR+both)/(LR+both+SR+Nocov))
    print('cov SR / covLR : ', (both)/(both+SR))
    print('(both+LR)/(both+SR) : ', (both+LR)/(both+SR))

    plt.figure(figsize=(4,4))
    v = venn3(subsets=(Nocov, 0, LR, 0, SR, 0, both), set_labels = ('Total Matrix', 'LR', 'SR'))
    v.get_patch_by_id('100').set_color('grey')
    v.get_patch_by_id('110').set_color('blue')
    v.get_patch_by_id('101').set_color('red')
    v.get_patch_by_id('111').set_color('purple')
    for text in v.set_labels:
        text.set_fontsize(15)
    for x in range(len(v.subset_labels)):
        if v.subset_labels[x] is not None:
            v.subset_labels[x].set_fontsize(15)
    plt.savefig('{}.mtrx_coverage_venn.png'.format(args.sample), dpi=600)

    ### Venn Tumor Only
    freq_SR = freq_SR[[i for i in tumor_cells if i in freq_SR.columns]]
    freq = freq[[i for i in tumor_cells if i in freq.columns]]

    venn = freq+freq_SR
    both = (venn==3).sum().sum()
    SR = (venn==2).sum().sum()
    LR = (venn==1).sum().sum()
    Nocov = (venn==0).sum().sum()

    print()
    print('Venn Diagram Values Cancer Only')
    print('Both :', both)
    print('LR only:', LR)
    print('SR only:', SR)
    print('Non covered:', Nocov)
    print('Perc cov LR : ', (LR+both)/(LR+both+SR+Nocov))
    print('Perc cov SR : ', (SR+both)/(LR+both+SR+Nocov))
    print('cov SR / covLR : ', (both)/(both+SR))
    print('(both+LR)/(both+SR) : ', (both+LR)/(both+SR))


    plt.figure(figsize=(4,4))
    v = venn3(subsets=(Nocov, 0, LR, 0, SR, 0, both), set_labels = ('Total Matrix', 'LR', 'SR'))
    v.get_patch_by_id('100').set_color('grey')
    v.get_patch_by_id('110').set_color('blue')
    v.get_patch_by_id('101').set_color('red')
    v.get_patch_by_id('111').set_color('purple')
    for text in v.set_labels:
        text.set_fontsize(15)
    for x in range(len(v.subset_labels)):
        if v.subset_labels[x] is not None:
            v.subset_labels[x].set_fontsize(15)
    plt.savefig('{}.mtrx_coverage_venn_tumor.png'.format(args.sample), dpi=600)

    ###scDNA
    cmap = plt.get_cmap('Reds', 100)
    cmap.set_bad('grey')
    muts['TumorColor'] = muts['scDNASupport_Tumor'].apply(lambda x: matplotlib.colors.rgb2hex(cmap(x)))
    muts['NonTumorColor'] = muts['scDNASupport_NonTumor'].apply(lambda x: matplotlib.colors.rgb2hex(cmap(x)))
    muts = muts[muts['INDEX'].isin(freq.index)]
    muts['INDEX'] = pd.Categorical(muts['INDEX'], list(set(freq.index)))
    muts = muts.sort_values(by = ['INDEX'])
    muts.to_csv('{}.BnpC_input.vcf'.format(args.sample), sep = '\t', index = False)
    scDNA = muts[['INDEX','scDNASupport_Tumor', 'scDNASupport_NonTumor','TumorColor', 'NonTumorColor','scDNA_Tumor', 'scDNA_NonTumor']]
    fusions=[]
    for f in fusion.index:
        fusions.append([f,0,0, 'violet','violet','0:0,0','0:0,0'])
    fusions = pd.DataFrame(fusions, columns = ['INDEX','scDNASupport_Tumor', 'scDNASupport_NonTumor','TumorColor', 'NonTumorColor','scDNA_Tumor','scDNA_NonTumor'])
    scDNA = pd.concat([scDNA,fusions])
    scDNA['INDEX'] = scDNA['INDEX'].map(VAF)
    scDNA.to_csv('{}.scDNA_support.tsv'.format(args.sample), sep = '\t', index = False)
    ###

    scDNAsup = muts[['scDNASupport_Tumor','scDNASupport_NonTumor']]
    scDNAsup = scDNAsup.dropna(axis = 0, how = 'all')
    scDNAsup = scDNAsup.replace(np.nan, 0)
    scDNAsup['Diff'] = scDNAsup['scDNASupport_Tumor'] - scDNAsup['scDNASupport_NonTumor']
    print('scDNA sup: ', len(scDNAsup[scDNAsup['Diff']>0.2])/len(scDNAsup), '%')

    end=time.time()
    print('Total time:',end-start) 

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
        '--dist_vcf', type=str, help='distal cells vcf.gz file'
    )
    parser.add_argument(
        '--mindepth', type=int, default=5,
        help='minimum depth in ouput'
    )
    parser.add_argument(
        '--minref', type=int, default=0,
        help='minimum reference allele reads'
    )
    parser.add_argument(
        '--minalt', type=int, default=3,
        help='minimum alternative allele reads'
    )
    parser.add_argument(
        '--mincells', type=int, default=8,
        help='minimum cancer cells with alternative allele reads'
    )
    parser.add_argument(
        '--maxmuts', type=int, default=300,
        help='max mutations for BnpC'
    )
    parser.add_argument(
        '--snv', action='store_true', help='snv only'
    )
    parser.add_argument(
        '--gnomAD', type=str, help='path to gnomAD SQL database, see https://pypi.org/project/gnomad-db/'
    )
    parser.add_argument(
        '--bam_tum_LR', type=str,
        help='Absolute or relative path(s) to tumor bam file'
    )
    parser.add_argument(
        '--bam_dist_LR', type=str,
        help='Absolute or relative path(s) to Nontumor bam file'
    )
    parser.add_argument(
        '--bam_tum_SR', type=str,
        help='Absolute or relative path(s) to tumor bam file'
    )
    parser.add_argument(
        '--bam_dist_SR', type=str,
        help='Absolute or relative path(s) to Nontumor bam file'
    )
    parser.add_argument(
        '--bam_scDNA', type=str,
        help='Absolute or relative path(s) to scDNA bam file'
    )
    parser.add_argument(
        '--ctypes', type=str,
        help='Absolute or relative path(s) to input barcode-to-ctype file'
    )
    parser.add_argument(
        '--fusions', type=str,
        help='Absolute or relative path(s) to fusions file'
    )
    parser.add_argument(
        '--fusions_SR', type=str,
        help='Absolute or relative path(s) to fusions file'
    )
    parser.add_argument(
        '--scDNA_clones', type=str,
        help='Absolute or relative path(s) to input scDNA clones file'
    )
    parser.add_argument(
        '--sample', type=str, default='',
        help='output csv file listing cells mutated'
    )
    parser.add_argument(
        '--blacklist', type=none_or_str, nargs='?', default=None,
        help='blacklist of loci'
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
