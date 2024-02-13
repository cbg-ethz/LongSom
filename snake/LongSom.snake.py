#use in conda env snakemake.32.7

workdir: '/cluster/work/bewi/members/dondia/projects/long_reads_tree/ctat_mut'
BINPATH = '/cluster/work/bewi/members/dondia/projects/long_reads_tree/bin'
REFPATH = '/cluster/work/bewi/members/dondia/projects/long_reads_tree/ref'
DATAPATH = '/cluster/work/bewi/members/dondia/projects/long_reads_tree/data'
MAXMUTS = config['filters']['maxmuts']
SUFFIX = config['BnpC']['suffix']
MISSDATA = config['BnpC']['missdata']

SAMPLES = ['B486_Tum','B486_Om','B497_Tum','B500_Tum','B500_Om']
PATIENTS = ['B486','B497','B500']

ruleorder: filter_vcf_with_dist > filter_vcf_no_dist

rule all:
    input:
        expand('{sample}_out/{sample}.filtered.vcf.gz', sample=SAMPLES),
        expand(f'{{patient}}_Tum_out/BnpC/{MAXMUTS}muts{SUFFIX}/genoCluster_posterior_mean.pdf', patient = PATIENTS),
        expand(f'{{patient}}_Tum_out/BnpC/{MAXMUTS}muts_SR{SUFFIX}/genoCluster_posterior_mean.pdf', patient = PATIENTS)


rule ctat_mutations:
    input:
        '../data/fastq/{sample}_rawdedup.test.fq',
    output:
        '{sample}_out/{sample}.filtered.vcf.gz'
    threads: 10
    resources:
        time = 1200,
        mem_mb=8000
    singularity:
        'ctat_mutations.v4.0.0.simg'
    shell:
        "/usr/local/src/ctat-mutations/ctat_mutations "
        "--left /data/fastq/{wildcards.sample}_rawdedup.test.fq "
        "--is_long_reads "
        "--output {wildcards.sample}_out "
        "--sample_id {wildcards.sample} "
        "--cpu 10 "
        "--genome_lib_dir /ref "

rule filter_vcf_with_dist:
    input:
        vcf_tum = '{patient}_Tum_out/{patient}_Tum.filtered.vcf.gz',
        vcf_dist = '{patient}_Om_out/{patient}_Om.filtered.vcf.gz',
        bam_tum = '{patient}_Tum_out/{patient}_Tum.bqsr.bam',
        bam_dist = '{patient}_Om_out/{patient}_Om.bqsr.bam'
    output:
        f'{{patient}}_Tum_out/{{patient}}_Tum.BnpC_input{MISSDATA}.{MAXMUTS}muts.tsv',
        f'{{patient}}_Tum_out/{{patient}}_Tum.BnpC_input_SR{MISSDATA}.{MAXMUTS}muts.tsv',
        f'{{patient}}_Tum_out/{{patient}}_Tum.scDNA_support.tsv'
    threads:8
    resources:
        time = 1200,
        mem_mb=8000
    params:
        bin_path=BINPATH,
        data_path=DATAPATH,
        ref_path=REFPATH,
        script = config['filters']['script'],
        maxmuts = MAXMUTS,
        minalt= config['filters']['minalt'],
        blacklist = lambda w: config['filters'][w.get("patient")]
    shell:
        "python {params.bin_path}/ctat_mut/create_input_bnpc_{params.script}.py "
        "--cpu {threads} --snv --tum_vcf {input.vcf_tum} "
        "--dist_vcf {input.vcf_dist} --minalt {params.minalt} --maxmuts {params.maxmuts} "
        "--sample {wildcards.patient}_Tum_out/{wildcards.patient}_Tum "
        "--ctypes ./ctypes/{wildcards.patient}_all.txt "
        "--fusions ./fusions/{wildcards.patient}.fusion_of_interest.tsv "
        "--fusions_SR ./fusions/{wildcards.patient}.fusion_of_interest.SR.tsv  "
        "--bam_tum_LR {input.bam_tum} --bam_dist_LR {input.bam_dist} "
        "--bam_tum_SR {params.data_path}/SR/{wildcards.patient}_Tum_SR.bam "
        "--bam_dist_SR {params.data_path}/SR/{wildcards.patient}_Om_SR.bam "
        "--bam_scDNA {params.data_path}/scDNA/{wildcards.patient}_scDNA.bam "
        "--scDNA_clones {params.data_path}/scDNA/clones_{wildcards.patient}.tsv "
        "--gnomAD {params.ref_path}/gnomAD --blacklist {params.blacklist}"

rule filter_vcf_no_dist:
    input:
        vcf_tum = '{patient}_Tum_out/{patient}_Tum.filtered.vcf.gz',
        bam_tum = '{patient}_Tum_out/{patient}_Tum.bqsr.bam'
    output:
        f'{{patient}}_Tum_out/{{patient}}_Tum.BnpC_input{MISSDATA}.{MAXMUTS}muts.tsv',
        f'{{patient}}_Tum_out/{{patient}}_Tum.BnpC_input_SR{MISSDATA}.{MAXMUTS}muts.tsv',
        f'{{patient}}_Tum_out/{{patient}}_Tum.scDNA_support.tsv'
    threads:8
    resources:
        time = 1200,
        mem_mb=8000
    params:
        bin_path = BINPATH,
        data_path = DATAPATH,
        ref_path=REFPATH,
        script = config['filters']['script'],
        maxmuts=MAXMUTS,
        minalt = config['filters']['minalt'],
        blacklist = lambda w: config['filters'][w.get("patient")]
    shell:
        "python {params.bin_path}/ctat_mut/create_input_bnpc_{params.script}.py "
        "--cpu {threads} --snv --tum_vcf {input.vcf_tum} "
        "--minalt {params.minalt} --maxmuts {params.maxmuts} "
        "--sample {wildcards.patient}_Tum_out/{wildcards.patient}_Tum "
        "--ctypes ./ctypes/{wildcards.patient}_all.txt "
        "--fusions ./fusions/{wildcards.patient}.fusion_of_interest.tsv  "
        "--fusions_SR ./fusions/{wildcards.patient}.fusion_of_interest.SR.tsv  "
        "--bam_tum_LR {input.bam_tum} "
        "--bam_tum_SR {params.data_path}/SR/{wildcards.patient}_Tum_SR.bam "
        "--bam_scDNA {params.data_path}/scDNA/{wildcards.patient}_scDNA.bam "
        "--scDNA_clones {params.data_path}/scDNA/clones_{wildcards.patient}.tsv "
        "--gnomAD {params.ref_path}/gnomAD --blacklist {params.blacklist}"

rule BnpC_clones:
    input:
        f'{{patient}}_Tum_out/{{patient}}_Tum.BnpC_input{MISSDATA}.{MAXMUTS}muts.tsv',
        f'{{patient}}_Tum_out/{{patient}}_Tum.scDNA_support.tsv'
    output:
        f'{{patient}}_Tum_out/BnpC/{MAXMUTS}muts{SUFFIX}/genoCluster_posterior_mean.pdf'
    threads:16
    resources:
        time = 1200,
        mem_mb=16000
    params:
        bin_path=BINPATH,
        data_path=DATAPATH,
        mcmc_steps = config['BnpC']['mcmc_steps'],
        maxmuts=MAXMUTS,
        suffix=SUFFIX,
        missdata=MISSDATA,
        dpa = lambda w: config['BnpC']['dpa'][w.get("patient")],
        cup = config['BnpC']['cup'],
        eup = config['BnpC']['eup'],
        FP = config['BnpC']['FP'],
        FN = config['BnpC']['FN'],
        estimator = config['BnpC']['estimator']
    conda:
        "BnpC"
    shell:
        "python {params.bin_path}/BnpC/run_BnpC.py "
        "{wildcards.patient}_Tum_out/{wildcards.patient}_Tum.BnpC_input{params.missdata}.{params.maxmuts}muts.tsv "
        "-n {threads} -s {params.mcmc_steps} -o {wildcards.patient}_Tum_out/BnpC/{params.maxmuts}muts{params.suffix} "
        "-e {params.estimator} -ap {params.dpa} -cup {params.cup} -FP {params.FP} -FN {params.FN} "
        "-pp 1 1 --ctypes {wildcards.patient}_Tum_out/{wildcards.patient}_Tum.celltype_correction.tsv "
        "--scDNA {wildcards.patient}_Tum_out/{wildcards.patient}_Tum.scDNA_support.tsv"

rule BnpC_clones_SR:
    input:
        f'{{patient}}_Tum_out/{{patient}}_Tum.BnpC_input_SR{MISSDATA}.{MAXMUTS}muts.tsv',
        f'{{patient}}_Tum_out/{{patient}}_Tum.scDNA_support.tsv'
    output:
        f'{{patient}}_Tum_out/BnpC/{MAXMUTS}muts_SR{SUFFIX}/genoCluster_posterior_mean.pdf'
    threads:16
    resources:
        time = 1200,
        mem_mb=16000
    params:
        bin_path=BINPATH,
        data_path=DATAPATH,
        mcmc_steps = config['BnpC']['mcmc_steps'],
        maxmuts=MAXMUTS,
        suffix=SUFFIX,
        missdata=MISSDATA,
        dpa = lambda w: config['BnpC']['dpa'][w.get("patient")],
        cup = config['BnpC']['cup'],
        eup = config['BnpC']['eup'],
        FP = config['BnpC']['FP'],
        FN = config['BnpC']['FN'],
        estimator = config['BnpC']['estimator']
    conda:
        "BnpC"
    shell:
        "python {params.bin_path}/BnpC/run_BnpC.py "
        "{wildcards.patient}_Tum_out/{wildcards.patient}_Tum.BnpC_input_SR{params.missdata}.{params.maxmuts}muts.tsv "
        "-n {threads} -s {params.mcmc_steps} -o {wildcards.patient}_Tum_out/BnpC/{params.maxmuts}muts_SR{params.suffix} "
        "-e {params.estimator} -ap {params.dpa} -cup {params.cup} -FP {params.FP} -FN {params.FN} "
        "-pp 1 1 --ctypes {wildcards.patient}_Tum_out/{wildcards.patient}_Tum.celltype_correction.tsv "
        "--scDNA {wildcards.patient}_Tum_out/{wildcards.patient}_Tum.scDNA_support.tsv"

