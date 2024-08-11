### This snakemake pipeline first finds the High Confidence Cancer variants (HCCV)
### Then cells are reannotated based on their SNV/fusion HCCV mutation status

import pandas as pd
CTYPES_REANNO = config['CellTypeReannotation']['celltypes']
OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']

include: 'SComaticPoN.smk'
include: 'CTATFusion.smk'

def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all:
    input:
        expand(f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv", id=IDS)
    default_target: True
    
rule Mapping:
    input:
        fastq = f'{DATA}/fastq/{{id}}.fastq.gz'
    output:
        sam = temp(f"{DATA}/bam/{{id}}.sam"),
    params:
        hg38 = config['Global']['genome']
    conda:
        "isoseq"
    threads:
        32
    resources:
        mem_mb = get_mem_mb
    shell:
        "minimap2 -t 30 -ax splice -uf --secondary=no -C5 "
        "{params.hg38} {input.fastq} > {output.sam}"

rule SortBam:
    input:
        sam = f"{DATA}/bam/{{id}}.sam"
    output:
        bam = temp(f"{DATA}/bam/{{id}}.bam"),
        bai = temp(f"{DATA}/bam/{{id}}.bam.bai")
    conda:
        "samtools"
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    shell:
        "samtools sort -@ {threads} {input.sam} -o {output.bam}##idx##{output.bai} --write-index"

rule AddBarcodeTag:
    input:
        bam = f"{DATA}/bam/{{id}}.bam",
        bai = f"{DATA}/bam/{{id}}.bam.bai"
    output:
        taggedbam = f"{DATA}/bam/{{id}}.CB.bam",
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python {params.scomatic}/AddBarcodeTag/AddBarcodeTag.py  "
        "--input {input.bam} --output {output.taggedbam} --cpu {threads}"

rule IndexBarcodedBam:
    input:
        taggedbam = f"{DATA}/bam/{{id}}.CB.bam",
    output:
        bai_bc = f"{DATA}/bam/{{id}}.CB.bam.bai"
    conda:
        "samtools"
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    shell:
        "samtools index -@ {threads}  {input.taggedbam}"

rule SplitBam_Reanno:
    input:
        bam = f"{DATA}/bam/{{id}}.CB.bam",
        bai = f"{DATA}/bam/{{id}}.CB.bam.bai",
        barcodes = f"{DATA}/ctypes/{{id}}.txt"
    output:
        expand("{OUTDIR}/CellTypeReannotation/SplitBam/{{id}}.{celltype}.bam", 
            celltype=CTYPES_REANNO, OUTDIR=[OUTDIR])
    resources:
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/SplitBam",
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality']
    shell:
        "python {params.scomatic}/SplitBam/SplitBamCellTypes.py  "
        "--bam {input.bam} --meta {input.barcodes} --id {wildcards.id} "
        "--outdir {params.outdir} --min_MQ {params.mapq}"

rule BaseCellCounter_Reanno:
    input:
        bam=f"{OUTDIR}/CellTypeReannotation/SplitBam/{{id}}.{{celltype}}.bam"
    output:
        tsv=f"{OUTDIR}/CellTypeReannotation/BaseCellCounter/{{id}}/{{id}}.{{celltype}}.tsv",
        tmp=temp(directory(f"{OUTDIR}/CellTypeReannotation/BaseCellCounter/{{id}}/temp_{{celltype}}/"))
    threads:
        32
    resources:
        time = 1200,
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        outdir=f"{OUTDIR}/CellTypeReannotation/BaseCellCounter",
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        chrom=config['SComatic']['BaseCellCounter']['chromosomes'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
    shell:
        "python {params.scomatic}/BaseCellCounter/BaseCellCounter.py "
        "--bam {input.bam} --ref {params.hg38} --chrom {params.chrom} "
        "--out_folder {params.outdir}/{wildcards.id}/ --nprocs {threads} "
        "--min_mq {params.mapq} --tmp_dir {output.tmp} "

rule MergeCounts_Reanno:
    input:
        expand("{OUTDIR}/CellTypeReannotation/BaseCellCounter/{{id}}/{{id}}.{celltype}.tsv", 
            celltype=CTYPES_REANNO, OUTDIR=[OUTDIR])
    output:
        tsv = f"{OUTDIR}/CellTypeReannotation/MergeCounts/{{id}}.BaseCellCounts.AllCellTypes.tsv"
    resources:
        time = 120,
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/BaseCellCounter",
    shell:
        "python {params.scomatic}/MergeCounts/MergeBaseCellCounts.py "
        "--tsv_folder {params.outdir}/{wildcards.id}/ --outfile {output.tsv}"

rule BaseCellCalling_step1_Reanno:
    input: 
        bb = f"{OUTDIR}/Pon/PoN/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/CellTypeReannotation/MergeCounts/{{id}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step1.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = get_mem_mb
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/BaseCellCalling",
        hg38=config['Global']['genome'],
        min_cell_types = config['SComatic']['BaseCellCalling']['Min_cell_types'],
        alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
        beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--ref  {params.hg38} --min_cell_types {params.min_cell_types} "
        "--alpha1 {params.alpha1} --beta1 {params.beta1} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} "

rule BaseCellCalling_step2_Reanno:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step1.tsv",
        pon_LR = f"{OUTDIR}/PoN/PoN/PoN_LR.tsv"
    output:
        f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step2.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/BaseCellCalling",
        RNA_editing = config['SComatic']['BaseCellCalling']['RNA_editing'],
        min_distance = config['SComatic']['BaseCellCalling']['min_distance'],
        pon_SR = config['SComatic']['BaseCellCalling']['PoN_SR'],
        gnomAD_db = config['SComatic']['BaseCellCalling']['gnomAD_db'],
        max_gnomAD_VAF = config['SComatic']['BaseCellCalling']['max_gnomAD_VAF'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step2.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--editing {params.RNA_editing} --min_distance {params.min_distance} "
        "--pon_SR {params.pon_SR} --pon_LR {input.pon_LR} "
        "--gnomAD_db {params.gnomAD_db} --gnomAD_max {params.max_gnomAD_VAF}"


rule HighConfidenceCancerVariants:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/CellTypeReannotation/HCCV/{{id}}.HCCV.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/HCCV",
        deltaVAF=config['SComatic']['HCCV']['deltaVAF'],
        deltaCCF=config['SComatic']['HCCV']['deltaCCF'],
        cancer = config['CellTypeReannotation']['cancer_ctype'],
    shell:
        "python {params.scomatic}/HighConfidenceCancerVariants/HighConfidenceCancerVariants.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer}"

rule HCCVSingleCellGenotype:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/HCCV/{{id}}.HCCV.tsv",
        bam = f"{DATA}/bam/{{id}}.CB.bam",
        barcodes = f"{DATA}/ctypes/{{id}}.txt",
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
    output:
        tsv=f"{OUTDIR}/CellTypeReannotation/HCCV/{{id}}.SingleCellGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/CellTypeReannotation/HCCV/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
        alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
        beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
        pval = config['SComatic']['SingleCellGenotype']['pvalue'],
        chrm_conta = config['SComatic']['chrM_contaminant']
    shell:
        "python {params.scomatic}/HighConfidenceCancerVariants/HCCVSingleCellGenotype.py "
        "--bam {input.bam} --infile {input.tsv} --outfile {output.tsv} "
        "--meta {input.barcodes} --alt_flag {params.alt_flag} --ref {params.hg38} "
        "--nprocs {threads} --min_mq {params.mapq} --pvalue {params.pval} "
        "--alpha1 {params.alpha1} --beta1 {params.beta1} "
        "--chrM_contaminant {params.chrm_conta} --tmp_dir {output.tmp}"

rule CellTypeReannotation:
    input:
        SNVs = f"{OUTDIR}/CellTypeReannotation/HCCV/{{id}}.SingleCellGenotype.tsv",
        fusions = f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv',
        barcodes = f"{DATA}/ctypes/{{id}}.txt"
    output:
        f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv"
    resources:
        time = 1200,
        mem_mb = get_mem_mb
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        min_variants = config['CellTypeReannotation']['min_variants']
    shell:
        "python {params.scomatic}/CellTypeReannotation/CellTypeReannotation.py "
        "--SNVs {input.SNVs} --fusions {input.fusions} --outfile {output} "
        "--meta {input.barcodes} --min_variants {params.min_variants} "


