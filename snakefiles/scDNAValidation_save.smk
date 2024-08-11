import pandas as pd 

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
SMPL=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']

all_clones_all_samples = {smpl : [] for smpl in SMPL}
input_betabin = {smpl : [] for smpl in SMPL}
input_basecellcounter= {smpl : [] for smpl in SMPL}

for smpl in SMPL:
    clones =  list(pd.read_csv(f"{DATA}/ctypes/scDNA/clones_{smpl}.tsv", sep='\t')['Cell_type'].unique())
    all_clones_all_samples[smpl]=clones
    for clone in clones:
        input_betabin[smpl].append(f"{OUTDIR}/scDNAValidation/BaseCellCounter/{smpl}/{smpl}.{clone}.tsv")
        input_basecellcounter[smpl].append(f"{OUTDIR}/scDNAValidation/BaseCellCounter/{smpl}/{smpl}.{clone}.tsv")
def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

def MergeCountsInput(wildcards):
    return expand("{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/{{scDNA}}.{clone}.tsv", 
            clone=all_clones_all_samples[wildcards.scDNA],
            OUTDIR=[OUTDIR])

def SplitBamOutput(wildcards):
    return expand("{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}.{clone}.bam", 
            clone=all_clones_all_samples[wildcards.scDNA],
            OUTDIR=[OUTDIR])

include: 'SNVCalling.smk'

rule all_scDNAValidation:
    input:
        f"{OUTDIR}/scDNAValidation/BetaBinEstimates.txt",
        expand(f"{OUTDIR}/scDNAValidation/BaseCellCalling/{{scDNA}}.calling.step1.tsv", scDNA = SMPL)
    default_target: True

rule SplitBam_scDNA:
    input:
        bam = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{scDNA}}.tsv"
    output:
        dynamic(f"{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}.{{clone}}.bam")
    resources:
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/SplitBam",
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality']
    shell:
        "python {params.scomatic}/SplitBam/SplitBamCellTypes.py  "
        "--bam {input.bam} --meta {input.barcodes} --id {wildcards.scDNA} "
        "--outdir {params.outdir} --min_MQ {params.mapq}"

rule BaseCellCounter_scDNA:
    input:
        bam = dynamic(f"{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}.{{clone}}.bam")
    output:
        tsv=f"{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/{{scDNA}}.{{clone}}.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/temp_{{clone}}/"))
    threads:
        32
    resources:
        time = 1200,
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        outdir=f"{OUTDIR}/scDNAValidation/BaseCellCounter",
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        chrom=config['SComatic']['BaseCellCounter']['chromosomes'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
        min_dp=config['scDNA']['BaseCellCounter']['min_dp'],
        min_cc=config['scDNA']['BaseCellCounter']['min_cc'],
    shell:
        "python {params.scomatic}/BaseCellCounter/BaseCellCounter.py "
        "--bam {input.bam} --ref {params.hg38} --chrom {params.chrom} "
        "--out_folder {params.outdir}/{wildcards.scDNA} "
        "--nprocs {threads} --tmp {output.tmp} --min_mq {params.mapq} "
        "--min_dp {params.min_dp} --min_cc {params.min_cc}"

rule CreateInputTsvListBetaBin_scDNA:
    input:
        input_betabin.values()
    output:
        temp(f"{OUTDIR}/scDNAValidation/BaseCellCounter/BaseCellCounter_files.txt")
    run:
        with open(output[0], "w") as out:
            for i in input:
                out.write(i+'\n')
    
rule BetaBinEstimation_scDNA:
    input:
        f"{OUTDIR}/scDNAValidation/BaseCellCounter/BaseCellCounter_files.txt"
    output:
        f"{OUTDIR}/scDNAValidation/BetaBinEstimates.txt"
    resources:
        time = 1200,
        mem_mb = 8000
    conda:
        "rpy2"
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python {params.scomatic}/BetaBinEstimation/BetaBinEstimation.py "
        "--in_tsv {input} --outfile {output}"

rule MergeCounts_scDNA:
    input:
        MergeCountsInput
    output:
        tsv = f"{OUTDIR}/scDNAValidation/MergeCounts/{{scDNA}}.BaseCellCounts.AllCellTypes.tsv"
    resources:
        time = 120,
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/BaseCellCounter",
    shell:
        "python {params.scomatic}/MergeCounts/MergeBaseCellCounts.py "
        "--tsv_folder {params.outdir}/{wildcards.scDNA} --outfile {output.tsv}"

rule BaseCellCalling_step1_scDNA:
    input: 
        bb = f"{OUTDIR}/scDNAValidation/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/scDNAValidation/MergeCounts/{{scDNA}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/scDNAValidation/BaseCellCalling/{{scDNA}}.calling.step1.tsv"
    resources:
        time = 120,
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/BaseCellCalling",
        hg38=config['Global']['genome'],
        alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
        beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        min_ac_reads = config['scDNA']['BaseCellCalling']['min_ac_reads'],
        min_ac_cells = config['scDNA']['BaseCellCalling']['min_ac_cells'],
        min_cells= config['scDNA']['BaseCellCalling']['min_cells'],
        min_cell_types = config['scDNA']['BaseCellCalling']['min_cell_types']
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.scDNA} "
        "--ref  {params.hg38} --alpha1 {params.alpha1} --beta1 {params.beta1} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} "
        "--min_ac_cells {params.min_ac_cells} --min_ac_reads {params.min_ac_reads} "
        "--min_cells {params.min_cells} --min_cell_types {params.min_cell_types}"


