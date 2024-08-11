import pandas as pd 

include: 'scDNACalling.smk'

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
SMPL=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']

all_clones_all_samples = {smpl : [] for smpl in SMPL}

for smpl in SMPL:
    clones =  list(pd.read_csv(f"{DATA}/ctypes/scDNA/clones_{smpl}.tsv", sep='\t')['Cell_type'].unique())
    all_clones_all_samples[smpl]=clones
print(all_clones_all_samples)

def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

def MergeCountsInput(wildcards):
    print(expand("{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/{{scDNA}}.{clone}.tsv", 
            clone=all_clones_all_samples[wildcards.scDNA],
            OUTDIR=[OUTDIR]))
    return expand("{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/{{scDNA}}.{clone}.tsv", 
            clone=all_clones_all_samples[wildcards.scDNA],
            OUTDIR=[OUTDIR])

def get_all_splitbam_file_names(wildcards):
    split_dir = str(checkpoints.SplitBam_scDNAValid.get(**wildcards).output["split"])
    CLONES = glob_wildcards(f'{split_dir}/{wildcards.scDNA}.{{clone}}.bam').clone
    print(CLONES)
    print(expand(f"{split_dir}/{wildcards.scDNA}.{{clone}}.bam",clone=CLONES))
    return expand(f"{split_dir}/{wildcards.scDNA}.{{clone}}.bam",clone=CLONES)



rule all_scDNAValidation:
    input:
        f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt",
        expand(f"{OUTDIR}/scDNAValidation/scDNAClonesGenotyping/{{scDNA}}.scDNAClonesGenotyping.tsv", 
        scDNA = SMPL),
        expand( f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step1.tsv",
        scDNA = SMPL)
    default_target: True

checkpoint SplitBam_scDNAValid:
    input:
        bam = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{scDNA}}.tsv"
    output:
        split = directory(f"{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}")
    resources:
        mem_mb = 1000
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/SplitBam",
        mapq=config['scDNA']['Validation']['BaseCellCounter']['min_mapping_quality']
    shell:
        "mkdir -p {output.split} && python {params.scomatic}/SplitBam/SplitBamCellTypes.py  "
        "--bam {input.bam} --meta {input.barcodes} --id {wildcards.scDNA} "
        "--outdir {params.outdir}/{wildcards.scDNA} --min_MQ {params.mapq}"

rule AggregateSplitBamOutput:
    input:
        bam=get_all_splitbam_file_names,
    output:
        touch(f"{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}.split_done.txt")

rule ExtractBed_scDNAValid:
    input:
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{scDNA}}.calling.step2.tsv"
    output:
        bed = f"{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}.scRNASites.bed"
    shell:
        """grep -v '#' {input.tsv} | cut -f1-3 > {output.bed}"""

rule BaseCellCounter_scDNAValid:
    input:
        txt=f"{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}.split_done.txt",
        bam=f"{OUTDIR}/scDNAValidation/SplitBam/{{scDNA}}/{{scDNA}}.{{clone}}.bam",
        bed=f"{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}.scRNASites.bed"
    output:
        tsv=f"{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/{{scDNA}}.{{clone}}.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/BaseCellCounter/{{scDNA}}/temp_{{clone}}/"))
    threads:
        32
    resources:
        time = 1200,
        mem_mb = 1000
    conda:
        "SComatic"
    params:
        outdir=f"{OUTDIR}/scDNAValidation/BaseCellCounter",
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        chrom=config['SComatic']['BaseCellCounter']['chromosomes'],
        mapq=config['scDNA']['Validation']['BaseCellCounter']['min_mapping_quality'],
        min_dp=config['scDNA']['Validation']['BaseCellCounter']['min_dp'],
        min_cc=config['scDNA']['Validation']['BaseCellCounter']['min_cc'],
    shell:
        "python {params.scomatic}/BaseCellCounter/BaseCellCounter.py "
        "--bam {input.bam} --ref {params.hg38} --chrom {params.chrom} "
        "--out_folder {params.outdir}/{wildcards.scDNA} --bed {input.bed} "
        "--nprocs {threads} --tmp {output.tmp} --min_mq {params.mapq} "
        "--min_dp {params.min_dp} --min_cc {params.min_cc}"

rule MergeCounts_scDNA:
    input:
        MergeCountsInput
    output:
        tsv = f"{OUTDIR}/scDNAValidation/MergeCounts/{{scDNA}}.BaseCellCounts.AllCellTypes.tsv"
    resources:
        time = 120,
        mem_mb = 4000
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/BaseCellCounter/Validation",
    shell:
        "python {params.scomatic}/MergeCounts/MergeBaseCellCounts.py "
        "--tsv_folder {params.outdir}/{wildcards.scDNA} --outfile {output.tsv}"

rule BaseCellCalling_step1_scDNA:
    input: 
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/scDNAValidation/MergeCounts/{{scDNA}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/scDNAValidation/BaseCellCalling/{{scDNA}}.calling.step1.tsv"
    resources:
        time = 120,
        mem_mb = 4000
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
        min_ac_reads = config['scDNA']['Validation']['BaseCellCalling']['min_ac_reads'],
        min_ac_cells = config['scDNA']['Validation']['BaseCellCalling']['min_ac_cells'],
        min_cells= config['scDNA']['Validation']['BaseCellCalling']['min_cells'],
        min_cell_types = config['scDNA']['Validation']['BaseCellCalling']['min_cell_types']
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.scDNA} "
        "--ref  {params.hg38} --alpha1 {params.alpha1} --beta1 {params.beta1} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} "
        "--min_ac_cells {params.min_ac_cells} --min_ac_reads {params.min_ac_reads} "
        "--min_cells {params.min_cells} --min_cell_types {params.min_cell_types}"

rule scDNAClonesGenotyping:
    input: 
        scDNA = f"{OUTDIR}/scDNAValidation/BaseCellCalling/{{scDNA}}.calling.step1.tsv",
        scRNA = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{scDNA}}.calling.step3.tsv"
    output:
        f"{OUTDIR}/scDNAValidation/scDNAClonesGenotyping/{{scDNA}}.scDNAClonesGenotyping.tsv"
    resources:
        time = 120,
        mem_mb = 1000
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/scDNAClonesGenotyping",
    shell:
        "python {params.scomatic}/scDNAClonesGenotyping/scDNAClonesGenotyping.py "
        "--scDNA {input.scDNA} --scDNA {input.scRNA} "
        "--outfile {params.outdir}/{wildcards.scDNA} "


