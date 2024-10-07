import pandas as pd 

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
NORMS=config['PoN']['ids']
SCOMATIC_PATH=config['Global']['scomatic']
NORM_CTYPES=config['PoN']['celltypes']

def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all_PoN:
    input:
        f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        f"{OUTDIR}/PoN/PoN/PoN.tsv"

rule mapping_PoN:
    input:
        fastq = f'{DATA}/fastq/{{norm}}.fastq.gz'
    output:
        sam = temp(f"{DATA}/bam/{{norm}}.sam"),
    wildcard_constraints:                                                                                                                                                            
        norm="([a-zA-Z]+)_Norm"   
    params:
        hg38 = config['Global']['genome']
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        mem_mb = get_mem_mb
    shell:
        "minimap2 -t 30 -ax splice -uf --secondary=no -C5 "
        "{params.hg38} {input.fastq} > {output.sam}"

rule sam_to_sortedbam_PoN:
    input:
        sam = f"{DATA}/bam/{{norm}}.sam"
    output:
        bam = temp(f"{DATA}/bam/{{norm}}.NoCB.bam"),
        bai = temp(f"{DATA}/bam/{{norm}}.NoCB.bam.bai")
    wildcard_constraints:                                                                                                                                                            
        norm="([a-zA-Z]+)_Norm"   
    conda:
        "envs/SComatic.yml"
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    shell:
        "samtools sort -@ {threads} {input.sam} -o {output.bam}##idx##{output.bai} --write-index"

rule AddBarcodeTag_PoN:
    input:
        bam = f"{DATA}/bam/{{norm}}.NoCB.bam",
        bai = f"{DATA}/bam/{{norm}}.NoCB.bam.bai"
    output:
        taggedbam = f"{DATA}/bam/{{norm}}.bam",
    wildcard_constraints:                                                                                                                                                            
        norm="([a-zA-Z]+)_Norm"    
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python {params.scomatic}/AddBarcodeTag/AddBarcodeTag.py  "
        "--input {input.bam} --output {output.taggedbam} --cpu {threads}"

rule IndexBarcodedBam_PoN:
    input:
        taggedbam = f"{DATA}/bam/{{norm}}.bam",
    output:
        bai_bc = f"{DATA}/bam/{{norm}}.bam.bai"
    wildcard_constraints:                                                                                                                                                            
        norm="([a-zA-Z]+)_Norm"     
    conda:
        "envs/SComatic.yml"
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    shell:
        "samtools index -@ {threads}  {input.taggedbam}"

rule SplitBam_PoN:
    input:
        bam = f"{DATA}/bam/{{norm}}.bam",
        bai = f"{DATA}/bam/{{norm}}.bam.bai",
        barcodes = f"{DATA}/ctypes/{{norm}}.txt"
    output:
        expand("{OUTDIR}/PoN/SplitBam/{{norm}}.{norm_celltype}.bam", 
            norm_celltype=NORM_CTYPES, OUTDIR=[OUTDIR])
    resources:
        mem_mb = 4096
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/PoN/SplitBam",
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality']
    shell:
        "python {params.scomatic}/SplitBam/SplitBamCellTypes.py  "
        "--bam {input.bam} --meta {input.barcodes} --id {wildcards.norm} "
        "--outdir {params.outdir} --min_MQ {params.mapq}"

rule BaseCellCounter_PoN:
    input:
        bam=f"{OUTDIR}/PoN/SplitBam/{{norm}}.{{norm_celltype}}.bam"
    output:
        tsv=f"{OUTDIR}/PoN/BaseCellCounter/{{norm}}/{{norm}}.{{norm_celltype}}.tsv",
        tmp=temp(directory(f"{OUTDIR}/PoN/BaseCellCounter/{{norm}}/temp_{{norm_celltype}}/"))
    threads:
        32
    resources:
        time = 1200,
        mem_mb = 1024
    conda:
        "envs/SComatic.yml"
    params:
        outdir=f"{OUTDIR}/PoN/BaseCellCounter",
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        chrom=config['SComatic']['BaseCellCounter']['chromosomes'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
    shell:
        "python {params.scomatic}/BaseCellCounter/BaseCellCounter.py "
        "--bam {input.bam} --ref {params.hg38} --chrom {params.chrom} "
        "--out_folder {params.outdir}/{wildcards.norm}/ "
        "--nprocs {threads} --min_mq {params.mapq} --tmp {output.tmp}"

rule CreateInputTsvListBetaBin_PoN:
    input:
        expand(f"{OUTDIR}/PoN/BaseCellCounter/{{norm}}/{{norm}}.{{norm_celltype}}.tsv", 
            norm=NORMS, norm_celltype=NORM_CTYPES)
    output:
        temp(f"{OUTDIR}/PoN/BaseCellCounter/BaseCellCounter_files.txt")
    run:
        with open(output[0], "w") as out:
            for i in input:
                out.write(i+'\n')
    
rule BetaBinEstimation_PoN:
    input:
        f"{OUTDIR}/PoN/BaseCellCounter/BaseCellCounter_files.txt"
    output:
        f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt"
    resources:
        time = 1200,
        mem_mb = 4096
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python {params.scomatic}/BetaBinEstimation/BetaBinEstimation.py "
        "--in_tsv {input} --outfile {output}"

rule MergeCounts_PoN:
    input:
        expand("{OUTDIR}/PoN/BaseCellCounter/{{norm}}/{{norm}}.{norm_celltype}.tsv", 
            norm_celltype=NORM_CTYPES, OUTDIR=[OUTDIR])
    output:
        tsv = f"{OUTDIR}/PoN/MergeCounts/{{norm}}.BaseCellCounts.AllCellTypes.tsv"
    resources:
        time = 120,
        mem_mb = 4096
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/PoN/BaseCellCounter",
    shell:
        "python {params.scomatic}/MergeCounts/MergeBaseCellCounts.py "
        "--tsv_folder {params.outdir}/{wildcards.norm}/ --outfile {output.tsv}"

rule BaseCellCalling_step1_PoN:
    input: 
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/PoN/MergeCounts/{{norm}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/PoN/BaseCellCalling/{{norm}}.calling.step1.tsv"
    resources:
        time = 120,
        mem_mb = 4096
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/PoN/BaseCellCalling",
        hg38=config['Global']['genome'],
        alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
        beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        min_ac_reads = config['PoN']['min_ac_reads'],
        min_ac_cells = config['PoN']['min_ac_cells'],
        min_cells= config['PoN']['min_cells'],
        min_cell_types = config['PoN']['min_cell_types']
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.norm} "
        "--ref  {params.hg38} --alpha1 {params.alpha1} --beta1 {params.beta1} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} "
        "--min_ac_cells {params.min_ac_cells} --min_ac_reads {params.min_ac_reads} "
        "--min_cells {params.min_cells} --min_cell_types {params.min_cell_types}"

rule CreateInputTsvPoN:
    input:
        expand(f"{OUTDIR}/PoN/BaseCellCalling/{{norm}}.calling.step1.tsv", 
            norm=NORMS)
    output:
        temp(f"{OUTDIR}/PoN/BaseCellCalling/BaseCellCalling_files.txt")
    run:
        with open(output[0], "w") as out:
            for i in input:
                out.write(i+'\n')

rule PoN:
    input:
        f"{OUTDIR}/PoN/BaseCellCalling/BaseCellCalling_files.txt"
    output:
        f"{OUTDIR}/PoN/PoN/PoN_LR.tsv"
    resources:
        time = 120,
        mem_mb = 4096
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python {params.scomatic}/PoN/PoN.py --in_tsv {input} "
        "--out_file {output} --min_samples 1 --rm_prefix No"
