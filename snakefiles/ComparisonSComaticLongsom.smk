### This snakemake pipeline detects SNVs before reannotation for comparison purposes
### It also checks what happens when unlocking all

import pandas as pd
CTYPES_REANNO = config['CellTypeReannotation']['celltypes']
OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']

#include: 'SNVCalling.smk'

def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all_ComparisonSComaticLongSom:
    input:
        expand(f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype/{{id}}.SingleCellGenotype.tsv",
         id=IDS),
        expand(f"{OUTDIR}/QC/BaseCellCalling/{{id}}.calling.step3.1000dist.tsv",
         id=IDS),
        expand(f"{OUTDIR}/QC/BaseCellCalling/{{id}}.calling.step3.NoPoNLR.tsv",
         id=IDS),
    default_target: True

rule BaseCellCalling_step3_Reanno:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['SNVCalling']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = config['SNVCalling']['clust_dist'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule SingleCellGenotype_Reanno:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/{{id}}.CB.bam",
        barcodes = f"{DATA}/ctypes/{{id}}.txt",
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        fusions = f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv',
    output:
        tsv=f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype/{{id}}.SingleCellGenotype.tsv",
        dp=f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype/{{id}}.DpMatrix.tsv",
        alt=f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype/{{id}}.AltMatrix.tsv",
        vaf=f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype/{{id}}.VAFMatrix.tsv",
        bin=f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype/{{id}}.BinaryMatrix.tsv",
        tmp=temp(directory(f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/SingleCellGenotype",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
        alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
        beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
        pval = config['SComatic']['SingleCellGenotype']['pvalue'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
    shell:
        "python {params.scomatic}/SingleCellGenotype/SingleCellGenotype.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} --fusions {input.fusions} "
        "--nprocs {threads} --min_mq {params.mapq} --pvalue {params.pval} "
        "--alpha1 {params.alpha1} --beta1 {params.beta1} --alt_flag {params.alt_flag} "
        "--chrM_contaminant {params.chrm_conta} --tmp_dir {output.tmp}"

rule BaseCellCalling_step3_SNVCalling_NoPoNLR:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/QC/BaseCellCalling/{{id}}.calling.step3.NoPoNLR.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/QC/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['SNVCalling']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = config['SNVCalling']['clust_dist'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py --PoN_LR False"
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule BaseCellCalling_step3_SNVCalling_1000dist:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/QC/BaseCellCalling/{{id}}.calling.step3.1000dist.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/QC/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['SNVCalling']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = 1000,
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule ComparisonSComaticLongSom:
    input:
        scomatic = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv",
        longsom = f"{OUTDIR}/SNVCallingn/BaseCellCalling/{{id}}.calling.step3.tsv"
    output:
        idk
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python params.scomatic}/ComparisonSComaticLongSom/ComparisonSComaticLongSom.py "
        ""
    

