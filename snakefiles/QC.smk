import pandas as pd

OUTDIR=config['Global']['outdir']
OUTDIR_SR = OUTDIR[:-2] + 'SR'
DATA=config['Global']['data']
IDS=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']
QC_PATH=config['Global']['qc']
CLUSTDIST = 100

include: 'SNVCalling.smk'

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all_QC:
    input:
        f"{OUTDIR}/QC/CompareNoPoN/PoNNoPoNComparison.png",
        f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/NoDistComparison.png",
        f"{OUTDIR}/QC/PlotCellTypeReannotation/MeanVAF.png",
        f"{OUTDIR}/QC/PlotCellTypeReannotation/MutationalBurden.png",
        f"{OUTDIR_SR}/Comparison/SRComparisonBoxplot.png",
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.UMAP.png", id = IDS),
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.Heatmap.png", id = IDS),
    default_target: True

rule Compare_PoN_NoPoN:
    input:
        LongSom = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        SComatic = f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step3.tsv",
        PoN = f"{OUTDIR}/PoN/PoN/PoN_LR.tsv"
    output:
        tsv = f"{OUTDIR}/QC/CompareNoPoN/{{id}}.CompareNoPoN.tsv"
    threads: 16
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/PoNComparison",
    shell:
        "python {params.qc}/PoNComparison/PoNComparison.py "
        "--LongSom {input.LongSom}  --SComatic {input.SComatic} "
        "--PoN {input.PoN} --id {wildcards.id} --outfile {output.tsv}"

rule plot_PoN_NoPoN:
    input:
        expand(f"{OUTDIR}/QC/CompareNoPoN/{{id}}.CompareNoPoN.tsv", id=IDS)
    output:
        png = f"{OUTDIR}/QC/CompareNoPoN/PoNNoPoNComparison.png"
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        indir=f"{OUTDIR}/QC/CompareNoPoN",
    shell:
        "python {params.qc}/PoNComparison/PlotPoNComparison.py "
        "--indir {params.indir} --outfile {output.png}"


rule BaseCellCalling_step3_QC_NoDist:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['SNVCalling']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = CLUSTDIST,
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule CloneGenotype_QC_NoDist:
    input: 
        tsv = f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv = f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/CloneGenotype/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/CloneGenotype/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/CloneGenotype",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
        beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        pval = config['SComatic']['SingleCellGenotype']['pvalue'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
    shell:
        "python {params.qc}/scDNAClonesGenotyping/scDNAClonesGenotyping.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} "
        "--alpha2 {params.alpha2} --beta2 {params.beta2} --pvalue {params.pval}  "
        "--chrM_contaminant {params.chrm_conta} --min_mq {params.mapq} "
        "--nprocs {threads} --tmp_dir {output.tmp}"


rule Compare_NoDist:
    input:
        geno = f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/CloneGenotype/{{id}}.CloneGenotype.tsv" if CLUSTDIST < 10000 else f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
        NoDist = f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/BaseCellCalling/{{id}}.calling.step3.tsv",
        YesDist = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv"
    output:
        tsv = f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/{{id}}.NoDistComparison.tsv"
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/",
    shell:
        "python {params.qc}/NoDistComparison/NoDistComparison.py "
        "--NoDist {input.NoDist}  --YesDist {input.YesDist} --id {wildcards.id} "
        "--geno {input.geno} --outfile {output.tsv}"

rule plot_NoDist:
    input:
        expand(f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/{{id}}.NoDistComparison.tsv", id=IDS)
    output:
        png = f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/NoDistComparison.png"
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        indir=f"{OUTDIR}/QC/NoDistComparison/{CLUSTDIST}/",
    shell:
        "python {params.qc}/NoDistComparison/PlotNoDistComparison.py "
        "--indir {params.indir} --outfile {output.png}"

rule PlotCellTypeReannotation:
    input:
        umap = f"{DATA}/UMAP/{{id}}.UMAP_coord.csv",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv",
        bin = f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.BinaryMatrix.tsv"
    output:
        umap = f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.UMAP.png",
        heatmap = f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.Heatmap.png",
        tsv = f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.boxplot.tsv"
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/PlotCellTypeReannotation",
    shell:
        "python {params.qc}/CellTypeReannotation/PlotCellTypeReannotation.py "
        "--umap {input.umap} --reannotation {input.barcodes} --binary {input.bin} "
        "--id {wildcards.id} --outfile {params.outdir}/{wildcards.id} "

rule PlotMutationalBurden:
    input:
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.boxplot.tsv", id=IDS)
    output:
        png1 = f"{OUTDIR}/QC/PlotCellTypeReannotation/MeanVAF.png",
        png2 = f"{OUTDIR}/QC/PlotCellTypeReannotation/MutationalBurden.png",
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        indir=f"{OUTDIR}/QC/PlotCellTypeReannotation",
    shell:
        "python {params.qc}/CellTypeReannotation/PlotMutationalBurden.py "
        "--indir {params.indir} --outfile1 {output.png1} --outfile2 {output.png2} "

rule CloneGenotype_QC_SR:
    input: 
        tsv = f"{OUTDIR_SR}/SComatic/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv = f"{OUTDIR_SR}/CellTypeReannotation/CloneGenotype/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR_SR}/CellTypeReannotation/CloneGenotype/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR_SR}/CellTypeReannotation/CloneGenotype/",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        pval = config['SComatic']['SingleCellGenotype']['pvalue'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
    shell:
        "python {params.qc}/scDNAClonesGenotyping/scDNAClonesGenotyping.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} "
        "--pvalue {params.pval}  "
        "--chrM_contaminant {params.chrm_conta} --min_mq {params.mapq} "
        "--nprocs {threads} --tmp_dir {output.tmp}"


rule Compare_LR_SR:
    input:
        geno_SR = f"{OUTDIR_SR}/CellTypeReannotation/CloneGenotype/{{id}}.CloneGenotype.tsv",
        geno_LR = f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
        SR = f"{OUTDIR_SR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv",
        LR = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv"
    output:
        tsv = f"{OUTDIR_SR}/Comparison/{{id}}.SRComparison.tsv"
    params:
        qc=QC_PATH,
    shell:
        "python {params.qc}/SRComparison/SRComparison.py "
        "--LR {input.LR}  --SR {input.SR} --id {wildcards.id} "
        "--geno_LR {input.geno_LR} --geno_SR {input.geno_SR} --outfile {output.tsv}"

rule plot_SRComparison:
    input:
        expand(f"{OUTDIR_SR}/Comparison/{{id}}.SRComparison.tsv", id=IDS)
    output:
        png = f"{OUTDIR_SR}/Comparison/SRComparisonBoxplot.png"
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        indir= f"{OUTDIR_SR}/Comparison/",
    shell:
        "python {params.qc}/SRComparison/PlotSRComparison.py "
        "--indir {params.indir} --outfile {output.png}"