#!/usr/bin/env python3

import argparse
import pandas as pd

def fusion_report(fusions,barcodes,min_ac_reads,min_ac_cells,max_MCF_noncancer,deltaMCF,outdir):
    barcodes = pd.read_table(barcodes)
    BC_cancer = barcodes[barcodes['Cell_type']=='Cancer']['Index']
    BC_noncancer = barcodes[barcodes['Cell_type']=='Non-Cancer']['Index']

    df = pd.read_table(fusions)
    df = df[df['SpliceType']=='ONLY_REF_SPLICE']
    df['List']=df['LR_accessions'].str.split(',')

    #Rename duplicates (different isoforms of the same gene-gene fusion)
    df['#FusionName'] = rename_duplicates(list(df['#FusionName']))

    # Exploding wide df into long df
    df_long = df.explode('List')
    df_long[['BC','UMI','ReadName']] = df_long['List'].str.split('^',expand=True)
    df_long = df_long[['#FusionName','LeftGene','LeftBreakpoint','RightGene','RightBreakpoint','SpliceType','BC','UMI','ReadName']]

    # num BC per fusion
    fusions=list(df['#FusionName'])

    n_bc_cancer={}
    n_bc_noncancer={}
    n_umi_cancer={}
    n_umi_noncancer={}
    for fusion in fusions:
        df_fus = df_long[df_long['#FusionName']==fusion]

        df_fus_cancer = df_fus[df_fus['BC'].isin(BC_cancer)]
        df_fus_noncancer = df_fus[df_fus['BC'].isin(BC_noncancer)]

        n_bc_cancer[fusion] = len(df_fus_cancer['BC'].unique())
        n_bc_noncancer[fusion] = len(df_fus_noncancer['BC'].unique())

        n_umi_cancer[fusion] = len(df_fus_cancer['UMI'].unique())
        n_umi_noncancer[fusion] = len(df_fus_noncancer['UMI'].unique())

    df['BC_Cancer'] = df['#FusionName'].map(n_bc_cancer)
    df['BC_Non-Cancer'] = df['#FusionName'].map(n_bc_noncancer)

    df['UMI_Cancer'] = df['#FusionName'].map(n_umi_cancer)
    df['UMI_Non-Cancer'] = df['#FusionName'].map(n_umi_noncancer)

    # Mutated cells fraction
    df['MCF_Cancer'] = df['BC_Cancer']/len(BC_cancer)
    df['MCF_Non-Cancer'] = df['BC_Non-Cancer']/len(BC_noncancer)

    df['Filter'] = df.apply(lambda x: filtering(x,min_ac_reads,min_ac_cells,max_MCF_noncancer,deltaMCF), axis=1)

    df = df[['#FusionName','Filter','UMI_Cancer','UMI_Non-Cancer','BC_Cancer','BC_Non-Cancer','MCF_Cancer','MCF_Non-Cancer','LeftGene','LeftLocalBreakpoint','LeftBreakpoint','RightGene','RightLocalBreakpoint','RightBreakpoint','SpliceType']]

    df.to_csv(outdir+'unfiltered.Fusions.tsv',sep='\t',index=False)

    df = df[df['Filter']=='PASS']

    df.to_csv(outdir+'.Fusions.tsv',sep='\t',index=False)

    df_long = df_long[df_long['#FusionName'].isin(df['#FusionName'])]
    df_long = df_long[['#FusionName','LeftGene','LeftBreakpoint','RightGene','RightBreakpoint','SpliceType','BC','UMI','ReadName']]

    df_long.to_csv(outdir+'.Fusions.SingleCellGenotype.tsv',sep='\t',index=False)


def filtering(x,min_ac_reads,min_ac_cells,max_MCF_noncancer,deltaMCF):
    if x['UMI_Cancer']<min_ac_reads:
        return 'Low_Cancer_UMI'
    if x['BC_Cancer']<min_ac_cells:
        return 'Low_Cancer_BC'
    if x['MCF_Non-Cancer']>0:
        if x['MCF_Cancer']-x['MCF_Non-Cancer']<deltaMCF:
            return 'Low_delta_MCF'
        if x['MCF_Non-Cancer']>max_MCF_noncancer:
            return 'High_Non-Cancer_MCF'
    return 'PASS'

def rename_duplicates(oldlist):
    newlist = []
    for i, v in enumerate(oldlist):
        totalcount = oldlist.count(v)
        count = oldlist[:i].count(v)
        newlist.append(v + str(count + 1) if totalcount > 1 else v)
    return newlist

def initialize_parser():
    parser = argparse.ArgumentParser(description='Report all fusion found in cancer cells and passing filters')
    parser.add_argument('--fusions', type=str, default=1, help='Fusions tsv (from ctat-fusion)', required = True)
    parser.add_argument('--barcodes', type=str, default=1, help='Barcodes to cell type (re)annotation tsv file', required = True)
    parser.add_argument('--min_ac_reads', type=int, default=3, help='Minimum fusion reads in cancer reads', required = True)
    parser.add_argument('--min_ac_cells', type=int, default=2, help='Minimum fusion reads in cancer cells', required = True)
    parser.add_argument('--max_MCF_noncancer', type=float, default=0.1, help='Maximum mutated cancer non-cancer cell fraction', required = True)
    parser.add_argument('--deltaMCF', type=float, default=0.3, help='Difference of mutated cell fraction between cancer and non-cancer cells', required = True)
    parser.add_argument('--outdir', type=str, default=1, help='Report tsv', required = True)
    return (parser)

def main():
    # 1. Arguments
    parser = initialize_parser()
    args = parser.parse_args()

    fusions=args.fusions
    barcodes=args.barcodes
    min_ac_reads=args.min_ac_reads
    min_ac_cells=args.min_ac_cells
    max_MCF_noncancer=args.max_MCF_noncancer
    deltaMCF=args.deltaMCF
    outdir=args.outdir

    # 2. Bam to fasta
    fusion_report(fusions,barcodes,min_ac_reads,min_ac_cells,max_MCF_noncancer,deltaMCF,outdir)


#-------------
# Execute code
#-------------

if __name__ == '__main__':
    main()
