### Fastq to barcoded bam

import pandas as pd
DATA=config['Global']['data']
IDS=config['Global']['ids']


rule all:
    input:
        expand(f"{DATA}/bam/{{id}}.bam.bai", id=IDS),
        expand(f"{DATA}/bam/{{id}}.bam", id=IDS),
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
        mem_mb = 1024
    shell:
        "minimap2 -t 30 -ax splice -uf --secondary=no -C5 "
        "{params.hg38} {input.fastq} > {output.sam}"

rule SortBam:
    input:
        sam = f"{DATA}/bam/{{id}}.sam"
    output:
        bam = temp(f"{DATA}/bam/{{id}}.NoCB.bam"),
        bai = temp(f"{DATA}/bam/{{id}}.NoCB.bam.bai")
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
        bam = f"{DATA}/bam/{{id}}.NoCB.bam",
        bai = f"{DATA}/bam/{{id}}.NoCB.bam.bai"
    output:
        taggedbam = f"{DATA}/bam/{{id}}.bam",
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
        taggedbam = f"{DATA}/bam/{{id}}.bam",
    output:
        bai_bc = f"{DATA}/bam/{{id}}.bam.bai"
    conda:
        "samtools"
    threads: 
        8
    resources:
        mem_mb = get_mem_mb
    shell:
        "samtools index -@ {threads}  {input.taggedbam}"