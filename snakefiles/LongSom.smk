include: 'BnpC.smk'
include: 'CellTypeReannotation.smk'
include: 'CTATFusion.smk'
#include: 'FastqBarcoding.smk'
include: 'InferCopyNumbers.smk'
include: 'SComaticPoN.smk'
include: 'SNVCalling.smk'
# Global
OUTDIR=config['Global']['outdir']
OUTBNPC=config['BnpC']['outdir']
IDS=config['Global']['ids']
# Run
REANNO=config['Run']['reanno']
SCOMATIC=config['Run']['scomatic']
ANNOVAR=config['Run']['annovar']
BNPC=config['Run']['bnpc']
CTATFUSION=config['Run']['ctatfusion']
INFERCNV=config['Run']['infercnv']
PLOT=config['Run']['plot']
SCDNA=config['Run']['scdna']

rule all:
    input:
        # Reannotated cell types
        expand(f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv", id=IDS) if REANNO else [],
        # Final SNV set
        expand(f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv", id=IDS) if SCOMATIC else [],
        # SNVs annotation 
        expand(f"{OUTDIR}/SNVCalling/Annotations/{{id}}.hg38_multianno.txt", id=IDS) if ANNOVAR else [],
        # Cell-Mutation matrix
        expand(f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.BinaryMatrix.tsv", id=IDS)if SCOMATIC else [],
        # Clonal reconstruction
        expand(f"{OUTDIR}/BnpC/{OUTBNPC}/{{id}}/genoCluster_posterior_mean_raw.pdf", id=IDS) if BNPC else [],
        # Fusions
        expand(f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv', id=IDS) if CTATFUSION else [],
        # CNVs
        expand(f"{OUTDIR}/InferCNV/InferCNV/{{id}}/infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.png", id = IDS) if INFERCNV else [],
        # Plots
        expand(f"{OUTDIR}/SNVCalling/ClusterMap/{{id}}.ClusterMap.Reannotation.pdf", id=IDS) if PLOT else [],
    default_target: True
