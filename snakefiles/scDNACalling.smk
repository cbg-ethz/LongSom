import pandas as pd 
include: 'SNVCalling.smk'

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
SMPL=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']

## scDNA specific
CLONE_TUMOR=config['scDNA']['Tumor_clone']
CLONE_NONTUMOR=config['scDNA']['NonTumor_clone']
CLONES = [CLONE_TUMOR,CLONE_NONTUMOR]

def get_mem_mb(wildcards, threads):
    return threads * 1024

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all_scDNACalling:
    input:
        f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt",
        expand(f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step3.tsv",
        scDNA = SMPL),
        expand(f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}.SingleCellGenotype.tsv",
        scDNA = SMPL),
    default_target: True

rule major_clone_only_scDNACalling:
    input:
        f"{DATA}/ctypes/scDNA/clones_{{scDNA}}.tsv"
    output:
        f"{DATA}/ctypes/scDNA/majorclones_{{scDNA}}.tsv"
    shell:
        "sed 's/Tum_[1-9]/Tum/g' {input} > {output}"

rule SplitBam_scDNACalling:
    input:
        bam = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/majorclones_{{scDNA}}.tsv"
    output:
        expand("{OUTDIR}/scDNACalling/SplitBam/{{scDNA}}.{clone}.bam", 
            clone=CLONES, OUTDIR=[OUTDIR])
    resources:
        mem_mb = 1000
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNACalling/SplitBam",
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality']
    shell:
        "python {params.scomatic}/SplitBam/SplitBamCellTypes.py  "
        "--bam {input.bam} --meta {input.barcodes} --id {wildcards.scDNA} "
        "--outdir {params.outdir} --min_MQ {params.mapq}"

# Selecting positions with enough coverage in scRNA
rule ExtractBedSNVCalling:
    input:
        tsv = ancient(f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{scDNA}}.calling.step1.tsv"),
    output:
        bed = f"{OUTDIR}/scDNACalling/BaseCellCounter/{{scDNA}}.scRNASites.bed"
    shell:
        """grep -v '#' {input.tsv} |"""
        """awk '{{if (($NF!="NA") && ($(NF-1)!="NA")) {{print $0}}}}' |"""
        """cut -f1-3 > {output.bed} """

rule BaseCellCounter_scDNACalling:
    input:
        bam=f"{OUTDIR}/scDNACalling/SplitBam/{{scDNA}}.{{clone}}.bam",
        bed = f"{OUTDIR}/scDNACalling/BaseCellCounter/{{scDNA}}.scRNASites.bed"
    output:
        tsv=f"{OUTDIR}/scDNACalling/BaseCellCounter/{{scDNA}}/{{scDNA}}.{{clone}}.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNACalling/BaseCellCounter/{{scDNA}}/temp_{{clone}}/"))
    threads:
        32
    resources:
        time = 1200,
        mem_mb = 1000
    conda:
        "SComatic"
    params:
        outdir=f"{OUTDIR}/scDNACalling/BaseCellCounter",
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        chrom=config['SComatic']['BaseCellCounter']['chromosomes'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        min_dp=config['scDNA']['BaseCellCounter']['min_dp'],
        min_cc=config['scDNA']['BaseCellCounter']['min_cc'],
    shell:
        "python {params.scomatic}/BaseCellCounter/BaseCellCounter.py "
        "--bam {input.bam} --ref {params.hg38} --chrom {params.chrom} "
        "--out_folder {params.outdir}/{wildcards.scDNA} --bed {input.bed} "
        "--nprocs {threads} --tmp_dir {output.tmp} --min_mq {params.mapq} "
        "--min_dp {params.min_dp} --min_cc {params.min_cc}"

rule CreateInputTsvListBetaBin_scDNACalling:
    input:
        expand(f"{OUTDIR}/scDNACalling/BaseCellCounter/{{scDNA}}/{{scDNA}}.{{clone}}.tsv", 
            scDNA = SMPL, clone=[CLONE_TUMOR])
    output:
        f"{OUTDIR}/scDNACalling/BaseCellCounter/BaseCellCounter_files.txt"
    run:
        with open(output[0], "w") as out:
            for i in input:
                out.write(i+'\n')
    
rule BetaBinEstimation_scDNACalling:
    input:
        f"{OUTDIR}/scDNACalling/BaseCellCounter/BaseCellCounter_files.txt"
    output:
        f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    resources:
        time = 1200,
        mem_mb = 10000
    conda:
        "rpy2"
    params:
        scomatic=SCOMATIC_PATH,
    shell:
        "python {params.scomatic}/BetaBinEstimation/BetaBinEstimation.py "
        "--in_tsv {input} --outfile {output}"

rule MergeCounts_scDNACalling:
    input:
        expand("{OUTDIR}/scDNACalling/BaseCellCounter/{{scDNA}}/{{scDNA}}.{clone}.tsv", 
            clone=CLONES, OUTDIR=[OUTDIR])
    output:
        tsv = f"{OUTDIR}/scDNACalling/MergeCounts/{{scDNA}}.BaseCellCounts.AllCellTypes.tsv"
    resources:
        time = 120,
        mem_mb = 4000
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNACalling/BaseCellCounter",
    shell:
        "python {params.scomatic}/MergeCounts/MergeBaseCellCounts.py "
        "--tsv_folder {params.outdir}/{wildcards.scDNA} --outfile {output.tsv}"

rule BaseCellCalling_step1_scDNACalling:
    input: 
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/scDNACalling/MergeCounts/{{scDNA}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step1.tsv"
    resources:
        time = 120,
        mem_mb = 4000
    conda:
        "SComatic"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNACalling/BaseCellCalling",
        hg38=config['Global']['genome'],
        alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
        beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        min_ac_reads = config['scDNA']['BaseCellCalling']['min_ac_reads'],
        min_ac_cells = config['scDNA']['BaseCellCalling']['min_ac_cells'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.scDNA} "
        "--ref  {params.hg38} --alpha1 {params.alpha1} --beta1 {params.beta1} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} "
        "--min_ac_reads {params.min_ac_reads} --min_ac_cells {params.min_ac_cells} "

rule BaseCellCalling_step2_scDNACalling:
    input: 
        tsv = f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step1.tsv",
        pon_LR = f"{OUTDIR}/PoN/PoN/PoN_LR.tsv"
    output:
        f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step2.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNACalling/BaseCellCalling",
        RNA_editing = config['SComatic']['BaseCellCalling']['RNA_editing'],
        min_distance = config['SComatic']['BaseCellCalling']['min_distance'],
        pon_SR = config['SComatic']['BaseCellCalling']['PoN_SR'],
        gnomAD_db = config['SComatic']['BaseCellCalling']['gnomAD_db'],
        max_gnomAD_VAF = config['SComatic']['BaseCellCalling']['max_gnomAD_VAF'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step2.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.scDNA} "
        "--editing {params.RNA_editing} --min_distance {params.min_distance} "
        "--pon_SR {params.pon_SR} --pon_LR {input.pon_LR} "
        "--gnomAD_db {params.gnomAD_db} --gnomAD_max {params.max_gnomAD_VAF}"


rule BaseCellCalling_step3_scDNACalling:
    input: 
        tsv = f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/scDNACalling/BaseCellCalling/{{scDNA}}.calling.step3.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNACalling/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['scDNA']['Tumor_clone'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['scDNA']['BaseCellCalling']['min_ac_reads'],
        clust_dist = config['SNVCalling']['clust_dist'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.scDNA} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule SingleCellGenotype_scDNACalling:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{scDNA}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{scDNA}}_scDNA.bam",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{scDNA}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv=f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}.SingleCellGenotype.tsv",
        dp=f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}.DpMatrix.tsv",
        alt=f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}.AltMatrix.tsv",
        vaf=f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}.VAFMatrix.tsv",
        bin=f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}.BinaryMatrix.tsv",
        tmp=directory(f"{OUTDIR}/scDNACalling/SingleCellGenotype/{{scDNA}}/")
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 1000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNACalling/SingleCellGenotype",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        pval = config['SComatic']['SingleCellGenotype']['pvalue'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
    shell:
        "python {params.scomatic}/SingleCellGenotype/SingleCellGenotypeNoFusion.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.scDNA} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} "
        "--nprocs {threads} --min_mq {params.mapq} --pvalue {params.pval} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} --alt_flag {params.alt_flag} "
        "--chrM_contaminant {params.chrm_conta} --tmp_dir {output.tmp}"
