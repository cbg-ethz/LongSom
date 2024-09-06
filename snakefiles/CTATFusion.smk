OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']


rule all_CTATFusion:
    input:
        expand(f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv', id=IDS)

    
rule CTATFusion:
    input:
        fastq = f'{DATA}/fastq/{{id}}.fastq.gz'
    output:
        f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv'
    threads: 32
    resources:
        time = 1200,
        mem_mb=97000
    singularity:
        'ctat_lr_fusion.v0.10.0.simg'
    shell:
        "ctat-LR-fusion "
        "-T /data/fastq/{wildcards.id}_rawdedup.test.fa "
        "--output {wildcards.id} "
        "--genome_lib_dir /ref "
        "--CPU {threads} --vis"
