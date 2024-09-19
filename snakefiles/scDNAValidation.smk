include: 'BnpC.smk'
include: 'QC.smk'
include: 'scDNACalling.smk'

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
SMPL=config['Global']['ids']
QC_PATH=config['Global']['qc']

def get_mem_mb(wildcards, threads):
    return threads * 1024

rule all_scDNAValidation:
    input:
        f"{OUTDIR}/scDNAValidation/Plots/F1_plot.png",
        f"{OUTDIR}/scDNAValidation/Plots/scDNAValidation.WaffleChart.LongSom.png",
        f"{OUTDIR}/scDNAValidation/Plots/scDNAValidation.WaffleChart.SComatic.png",
        expand(f"{OUTDIR}/scDNAValidation/Plots/{{id}}.Venn2.png",
         id = SMPL),
        expand(f"{OUTDIR}/scDNAValidation/Plots/{{id}}.Venn3.png",
         id = SMPL),
        #SNVCalling
        expand(f"{OUTDIR}/SNVCalling/ClusterMap/{{id}}.ClusterMap.Reannotation.pdf",
         id=IDS),
        # BnpC
        expand(f"{OUTDIR}/BnpC/{OUTBNPC}/{{id}}/genoCluster_posterior_mean_raw.pdf",
         id=IDS),
        # Reanno 
         expand(f"{OUTDIR}/SComatic/SingleCellGenotype/{{id}}.BinaryMatrix.tsv", 
         id=IDS),
        # QC
        f"{OUTDIR}/QC/CompareNoPoN/PoNNoPoNComparison.png",
        f"{OUTDIR}/QC/PlotCellTypeReannotation/MeanVAF.png",
        f"{OUTDIR}/QC/PlotCellTypeReannotation/MutationalBurden.png",
        f"{OUTDIR_SR}/Comparison/SRComparisonBoxplot.png",
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.UMAP.png", id = IDS),
        expand(f"{OUTDIR}/QC/PlotCellTypeReannotation/{{id}}.Heatmap.png", id = IDS),
    default_target: True

rule CloneGenotype_LongSom:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv = f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom",
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

rule CloneGenotype_SComatic:
    input: 
        tsv = f"{OUTDIR}/SComatic/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv",
        bb = f"{OUTDIR}/scDNACalling/BetaBinEstimates.txt"
    output:
        tsv = f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic",
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

rule scDNA_supp_in_scRNA:
    input: 
        tsv = f"{OUTDIR}/scDNACalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/{{id}}.bam",
        bai = f"{DATA}/bam/{{id}}.bam.bai",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv",
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
    output:
        tsv = f"{OUTDIR}/scDNAValidation/CloneGenotype/scDNASupportInRNA/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/CloneGenotype/scDNASupportInRNA/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/CloneGenotype/scDNASupportInRNA",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        ctype=config['CellTypeReannotation']['ctype_column'],
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
        "--nprocs {threads} --tmp_dir {output.tmp} --ctype {params.ctype}"

rule ComparisonSComaticLongSom:
    input:
        scDNACalls = f"{OUTDIR}/scDNACalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        scDNAValidLong = f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
        scDNAValidSCom = f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}.CloneGenotype.tsv",
        scDNA_supp_in_scRNA = f"{OUTDIR}/scDNAValidation/CloneGenotype/scDNASupportInRNA/{{id}}.CloneGenotype.tsv",
    output:
        f"{OUTDIR}/scDNAValidation/Plots/{{id}}.Venn2.png",
        f"{OUTDIR}/scDNAValidation/Plots/{{id}}.Venn3.png",
        f"{OUTDIR}/scDNAValidation/Plots/{{id}}.F1Scores.tsv",
    conda:
        "envs/BnpC.yml"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/Plots",
        dp=config['scDNAValidation']['min_depth'],
    shell:
        "python {params.qc}/ComparisonSComaticLongSom/ComparisonSComaticLongSom.py "
        "--id {wildcards.id}  --scDNACalls {input.scDNACalls} --dp {params.dp} "
        "--scDNAValidLong {input.scDNAValidLong} --scDNAValidSCom {input.scDNAValidSCom} "
        "--scDNA_supp_in_scRNA {input.scDNA_supp_in_scRNA} --outfile {params.outdir}/{wildcards.id} "

rule plot_ComparisonSComaticLongSom:
    input:
        expand(f"{OUTDIR}/scDNAValidation/Plots/{{id}}.F1Scores.tsv",
         id = SMPL),
    output:
        png1 = f"{OUTDIR}/scDNAValidation/Plots/F1_plot.png",
        png2 = f"{OUTDIR}/scDNAValidation/Plots/scDNASupport.png"
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        indir=f"{OUTDIR}/scDNAValidation/Plots",
    shell:
        "python {params.qc}/ComparisonSComaticLongSom/PlotComparisonSComaticLongSom.py "
        "--indir {params.indir} "

rule plot_scDNAValidation:
    input:
        expand(f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
         id = SMPL),
        expand(f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}.CloneGenotype.tsv",
         id = SMPL),
    output:
        png1 = f"{OUTDIR}/scDNAValidation/Plots/scDNAValidation.WaffleChart.LongSom.png",
        png2 = f"{OUTDIR}/scDNAValidation/Plots/scDNAValidation.WaffleChart.SComatic.png"
    conda:
        "envs/BnpC.yml"
    params:
        qc=QC_PATH,
        indir1=f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom",
        indir2=f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic",
        dp=config['scDNAValidation']['min_depth'],
    shell:
        "python {params.qc}/scDNAValidation/PlotWaffleChartscDNAValid.py "
        "--indir1 {params.indir1} --indir2 {params.indir2} "
        "--outfile1 {output.png1} --outfile2 {output.png2} "
        "--dp {params.dp}"

    

