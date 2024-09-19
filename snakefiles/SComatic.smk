### This snakemake pipeline applies OG SComatic

import pandas as pd
OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']
SCOMATIC_PATH= '/cluster/work/bewi/members/dondia/projects/long_reads_tree/bin/SComatic/scripts'
CTATFUSION=config['Run']['ctatfusion']

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all_reanno:
    input:
         expand(f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step3.tsv", id=IDS),
         expand(f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.SingleCellGenotype.tsv", id=IDS),

rule BaseCellCalling_step1_SComatic:
    input: 
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/CellTypeReannotation/MergeCounts/{{id}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step1.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 4096
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SComatic/BaseCellCalling",
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

rule BaseCellCalling_step2_SComatic:
    input: 
        tsv = f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step1.tsv",
        pon_LR = f"{OUTDIR}/PoN/PoN/{{id}}.PoN_LR.tsv"
    output:
        f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step2.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 4096
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SComatic/BaseCellCalling",
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

rule BaseCellCalling_step3_SComatic:
    input: 
        tsv = f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 4096
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SComatic/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['CellTypeReannotation']['cancer_ctype'],
        chrm_conta = config['CellTypeReannotation']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = config['SNVCalling']['clust_dist'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--clust_dist {params.clust_dist} "

rule SingleCellGenotype_SComatic:
    input: 
        tsv = f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/{{id}}.bam",
        barcodes = f"{DATA}/ctypes/{{id}}.txt",
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        fusions = f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv' if CTATFUSION else [],
    output:
        tsv=f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.SingleCellGenotype.tsv",
        dp=f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.DpMatrix.tsv",
        alt=f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.AltMatrix.tsv",
        vaf=f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.VAFMatrix.tsv",
        bin=f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.BinaryMatrix.tsv",
        tmp=temp(directory(f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 4096
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SComatic/SingleCellGenotype",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        pval = config['SComatic']['SingleCellGenotype']['pvalue'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
    shell:
        "python {params.scomatic}/SingleCellGenotype/SingleCellGenotype.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} --fusions {input.fusions} "
        "--nprocs {threads} --min_mq {params.mapq} --pvalue {params.pval} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} --alt_flag {params.alt_flag} "
        "--chrM_contaminant {params.chrm_conta} --tmp_dir {output.tmp}"
