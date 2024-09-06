import pandas as pd
CTYPES = config['SNVCalling']['celltypes']
OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']
CTATFUSION =config['Run']['ctatfusion']

#include: 'CellTypeReannotation.smk'

rule all_SNVCalling:
    input:
        expand(f"{OUTDIR}/SNVCalling/ClusterMap/{{id}}.ClusterMap.Reannotation.pdf",
         id=IDS),
        expand(f"{OUTDIR}/SNVCalling/Annotations/{{id}}.hg38_multianno.txt", 
         id=IDS),

rule SplitBam:
    input:
        bam = f"{DATA}/bam/{{id}}.bam",
        bai = f"{DATA}/bam/{{id}}.bam.bai",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv",
    output:
        expand("{OUTDIR}/SNVCalling/SplitBam/{{id}}.{celltype}.bam", 
            celltype=CTYPES, OUTDIR=[OUTDIR])
    resources:
        mem_mb = 32000
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/SplitBam",
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality']
    shell:
        "python {params.scomatic}/SplitBam/SplitBamCellTypes.py  "
        "--bam {input.bam} --meta {input.barcodes} --id {wildcards.id} "
        "--outdir {params.outdir} --min_MQ {params.mapq}"

rule BaseCellCounter:
    input:
        bam=f"{OUTDIR}/SNVCalling/SplitBam/{{id}}.{{celltype}}.bam"
    output:
        tsv=f"{OUTDIR}/SNVCalling/BaseCellCounter/{{id}}/{{id}}.{{celltype}}.tsv",
        tmp=temp(directory(f"{OUTDIR}/SNVCalling/BaseCellCounter/{{id}}/temp_{{celltype}}/"))
    threads:
        32
    resources:
        time = 1200,
        mem_mb = 1000
    conda:
        "envs/SComatic.yml"
    params:
        outdir=f"{OUTDIR}/SNVCalling/BaseCellCounter",
        scomatic=SCOMATIC_PATH,
        hg38=config['Global']['genome'],
        chrom=config['SComatic']['BaseCellCounter']['chromosomes'],
        mapq=config['SComatic']['BaseCellCounter']['min_mapping_quality'],
    shell:
        "python {params.scomatic}/BaseCellCounter/BaseCellCounter.py "
        "--bam {input.bam} --ref {params.hg38} --chrom {params.chrom} "
        "--out_folder {params.outdir}/{wildcards.id}/ --nprocs {threads} "
        "--min_mq {params.mapq} --tmp_dir {output.tmp} "

rule MergeCounts:
    input:
        expand("{OUTDIR}/SNVCalling/BaseCellCounter/{{id}}/{{id}}.{celltype}.tsv", 
            celltype=CTYPES, OUTDIR=[OUTDIR])
    output:
        tsv = f"{OUTDIR}/SNVCalling/MergeCounts/{{id}}.BaseCellCounts.AllCellTypes.tsv"
    resources:
        time = 120,
        mem_mb = 8000
    conda:
        "envs/SComatic.yml"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/BaseCellCounter",
    shell:
        "python {params.scomatic}/MergeCounts/MergeBaseCellCounts.py "
        "--tsv_folder {params.outdir}/{wildcards.id}/ --outfile {output.tsv}"

rule BaseCellCalling_step1:
    input: 
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        tsv = f"{OUTDIR}/SNVCalling/MergeCounts/{{id}}.BaseCellCounts.AllCellTypes.tsv"
    output:
        f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step1.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/BaseCellCalling",
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

rule BaseCellCalling_step2:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step1.tsv",
        pon_LR = f"{OUTDIR}/PoN/PoN/PoN_LR.tsv"
    output:
        f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step2.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/BaseCellCalling",
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


rule BaseCellCalling_step3:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "envs/SComatic.yml"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/BaseCellCalling",
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

rule SingleCellGenotype:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/{{id}}.bam",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv",
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
        fusions = f'{OUTDIR}/CTATFusion/{{id}}.fusion_of_interest.tsv' if CTATFUSION else [],
    output:
        tsv=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.SingleCellGenotype.tsv",
        dp=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.DpMatrix.tsv",
        alt=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.AltMatrix.tsv",
        vaf=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.VAFMatrix.tsv",
        bin=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.BinaryMatrix.tsv",
        tmp=temp(directory(f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 1000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/SingleCellGenotype",
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

rule clustermap:
    input:
        bin = f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.BinaryMatrix.tsv",
        ctypes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv"
    output:
        f"{OUTDIR}/SNVCalling/ClusterMap/{{id}}.ClusterMap.Reannotation.pdf",
        f"{OUTDIR}/SNVCalling/ClusterMap/{{id}}.ClusterMap.NoReannotation.pdf",
    conda:
        "envs/BnpC.yml"
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/SNVCalling/ClusterMap",
        height = config['SNVCalling']['clustermap']['height'],
        width = config['SNVCalling']['clustermap']['width'],
    shell:
        "python {params.scomatic}/ClusterMap/ClusterMap.py "
        "--bin {input.bin} --ctypes {input.ctypes} "
        "--height {params.height} --width {params.width} "
        "--out_dir {params.outdir}/{wildcards.id}"

rule input_annovar:
    input:
        f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
    output:
        f"{OUTDIR}/SNVCalling/Annotations/{{id}}.avinput" 
    shell:
        """grep -v '#' {input} |  tr '\t' '-' | """
        """awk -F'-' -v OFS='\t' '{{print $1,$2,$3,$4,$5,$0}}' > {output} """

rule run_annovar:
    input:
        f"{OUTDIR}/SNVCalling/Annotations/{{id}}.avinput"
    output:
        f"{OUTDIR}/SNVCalling/Annotations/{{id}}.hg38_multianno.txt"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        annovar = config['Global']['annovar'],
        outdir = f"{OUTDIR}/SNVCalling/Annotations/",
    shell:
        "perl {params.annovar}/table_annovar.pl {input} "
        "-out {params.outdir}/{wildcards.id} {params.annovar}/humandb/ "
        "-buildver hg38 -remove -nastring . -polish "
        "-xref {params.annovar}/example/gene_xref.txt "
        "-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a "
        "-operation gx,r,f,f,f "
