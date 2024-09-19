import pandas as pd

workdir: config['specific']['workdir']
SAMPLES = config['samples'].keys()
BIN = config['specific']['scripts']

def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_cellranger_rawfiles(wildcards):
    return '{}/raw_feature_bc_matrix/barcodes.tsv.gz'.format(config['samples'][wildcards.sample]["cellranger_folder"])

def get_cellranger_filteredfiles(wildcards):
    return '{}/filtered_feature_bc_matrix/barcodes.tsv'.format(config['samples'][wildcards.sample]["cellranger_folder"])

def sample2ids(wildcards):
    return expand('input_flntc/{{sample}}_{id}.fltnc.bam', 
                id = config['samples'][wildcards.sample]['ids'])    

rule all:
    input:
        expand('results/{sample}/fraction_alt_emptydrops.done',sample=SAMPLES)
        

rule identify_emptydroplets:
    input:
        raw = get_cellranger_rawfiles,
        filtered = get_cellranger_filteredfiles
    output:
        emptydrops = "emptydroplets/barcodes/{sample}.tsv"
    params:
        bin_path = BIN
    shell:
        "python {params.bin_path}/empty_droplets_listing.py --raw {input.raw} "
        "--filtered {input.filtered} --sample {wildcards.sample}"


rule emptydroplets_fastq:
    input:
        spl = sample2ids,
        emptydrops = "emptydroplets/barcodes/{sample}.tsv" 
    output:
        fastq = 'emptydroplets/{sample}.fastq.gz'
    params:
        bin_path = BIN
    conda:
        "pysam"
    threads:
        8
    resources:
        time = 1200,
        mem_mb = 32000
    shell:
        "python {params.bin_path}/empty_droplets_fastq.py --sample {wildcards.sample} "
        "--emptydroplets {input.emptydrops} --cpu {threads} "


rule mapping:
    input:
        fastq = 'emptydroplets/{sample}.fastq.gz'
    output:
        sam = "emptydroplets/{sample}.sam",
    params:
        hg38 = config['specific']['genome']
    conda:
        "isoseq"
    threads:
        32
    resources:
        mem_mb = get_mem_mb
    shell:
        "minimap2 -t 30 -ax splice -uf --secondary=no -C5 "
        "{params.hg38} {input.fastq} > {output.sam}"

rule sam_to_sortedbam:
    input:
        sam = ancient("emptydroplets/{sample}.sam")
    output:
        bam = "emptydroplets/{sample}.bam",
        bai = "emptydroplets/{sample}.bam.bai"
    conda:
        "samtools"
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    shell:
        "samtools sort -@ {threads} {input.sam} -o {output.bam}##idx##{output.bai} --write-index"

rule fraction_alt_emptydrops:
    input:
        bam = "emptydroplets/{sample}.bam",
        bai = "emptydroplets/{sample}.bam.bai",
        vcf = "longsom_muts/{sample}.BnpC_input.vcf"
    output:
        touch('results/{sample}/fraction_alt_emptydrops.done')
    conda:
        "pysam"
    threads:
        8
    resources:
        time = 1200,
        mem_mb = 32000
    params:
        bin_path = BIN
    shell:
        "python {params.bin_path}/fraction_alt_emptydrops.py --bam {input.bam} "
        "--vcf {input.vcf} --sample {wildcards.sample} --cpu {threads}"


