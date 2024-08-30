import pandas as pd

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']
QC_PATH=config['Global']['qc']

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all_QC:
    input:
        f"{OUTDIR}/QC/CompareNoPoN/PoNNoPoNComparison.png",
        f"{OUTDIR}/QC/PlotCellTypeReannotation/MutationalBurden.png",
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.UMAP.png", id = IDS),
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.Heatmap.png", id = IDS)
    default_target: True

rule ModifyBaseCellCalling_step2_QC_NoPoNLR:
    input:
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        tsv = f"{OUTDIR}/QC/CompareNoPoN/BaseCellCalling/{{id}}.calling.step2.tsv"
    shell:
        "sed 's/PoN_LR/PASS/g' {input.tsv} > {output.tsv}"

rule BaseCellCalling_step3_QC_NoPoNLR:
    input: 
        tsv = f"{OUTDIR}/QC/CompareNoPoN/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/QC/CompareNoPoN/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/QC/CompareNoPoN/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['SNVCalling']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = config['SNVCalling']['clust_dist'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule CloneGenotype_QC_NoPoNLR:
    input: 
        tsv = f"{OUTDIR}/QC/CompareNoPoN/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv = f"{OUTDIR}/QC/CompareNoPoN/CloneGenotype/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/QC/CompareNoPoN/CloneGenotype/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/CompareNoPoN/CloneGenotype",
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

rule Compare_PoN_NoPoN:
    input:
        geno = f"{OUTDIR}/QC/CompareNoPoN/CloneGenotype/{{id}}.CloneGenotype.tsv",
        NoPoN = f"{OUTDIR}/QC/CompareNoPoN/BaseCellCalling/{{id}}.calling.step3.tsv",
        YesPoN = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv"
    output:
        tsv = f"{OUTDIR}/QC/CompareNoPoN/{{id}}.CompareNoPoN.tsv"
    conda:
        "BnpC"
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/PoNComparison",
    shell:
        "python {params.qc}/PoNComparison/PoNComparison.py "
        "--NoPoN {input.NoPoN}  --YesPoN {input.YesPoN} --id {wildcards.id} "
        "--genotype {input.geno} --outfile {output.tsv}"

rule plot_PoN_NoPoN:
    input:
        expand(f"{OUTDIR}/QC/CompareNoPoN/{{id}}.CompareNoPoN.tsv", id=IDS)
    output:
        png = f"{OUTDIR}/QC/CompareNoPoN/PoNNoPoNComparison.png"
    conda:
        "BnpC"
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
        f"{OUTDIR}/QC/NoDistComparison/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/QC/NoDistComparison/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['SNVCalling']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = 0,
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule CloneGenotype_QC_NoDist:
    input: 
        tsv = f"{OUTDIR}/QC/NoDistComparison/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv = f"{OUTDIR}/QC/NoDistComparison/CloneGenotype/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/QC/NoDistComparison/CloneGenotype/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/NoDistComparison/CloneGenotype",
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
        geno = f"{OUTDIR}/QC/NoDistComparison/CloneGenotype/{{id}}.CloneGenotype.tsv",
        NoDist = f"{OUTDIR}/QC/NoDistComparison/BaseCellCalling/{{id}}.calling.step3.tsv",
        YesDist = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv"
    output:
        tsv = f"{OUTDIR}/QC/NoDistComparison/{{id}}.NoDistComparison.tsv"
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/QC/NoDistComparison",
    shell:
        "python {params.qc}/NoDistComparison/NoDistComparison.py "
        "--NoDist {input.NoPoN}  --YesDist {input.YesPoN} --id {wildcards.id}"
        "--genotype {input.geno} --outfile {output.tsv}"

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
        "BnpC"
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
        png = f"{OUTDIR}/QC/PlotCellTypeReannotation/MutationalBurden.png",
    conda:
        "BnpC"
    params:
        qc=QC_PATH,
        indir=f"{OUTDIR}/QC/PlotCellTypeReannotation",
    shell:
        "python {params.qc}/CellTypeReannotation/PlotMutationalBurden.py "
        "--indir {params.indir} --outfile {output.png} "