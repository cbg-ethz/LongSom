include: 'scDNACalling.smk'

OUTDIR=config['Global']['outdir']
DATA=config['Global']['data']
SMPL=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']

def get_mem_mb(wildcards, threads):
    return threads * 1024

rule all_scDNAValidation:
    input:
        expand(f"{OUTDIR}/scDNAValidation/Venn3/{{id}}.Venn3.png",
        id = SMPL),
        expand(f"{OUTDIR}/scDNAValidation/Venn3/{{id}}.LongSomScores.tsv",
        id = SMPL),
    default_target: True

rule BaseCellCalling_step3_Reanno:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step2.tsv"
    output:
        f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv"
    conda:
        "SComatic"
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/CellTypeReannotation/BaseCellCalling",
        deltaVAF=config['SComatic']['BaseCellCalling']['deltaVAF'],
        deltaCCF=config['SComatic']['BaseCellCalling']['deltaCCF'],
        cancer = config['CellTypeReannotation']['cancer_ctype'],
        chrm_conta = config['SComatic']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['min_ac_reads'],
        clust_dist = config['SNVCalling']['clust_dist'],
    shell:
        "python {params.scomatic}/BaseCellCalling/BaseCellCalling.step3.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} --chrM_contaminant {params.chrm_conta} "
        "--deltaVAF {params.deltaVAF} --deltaCCF {params.deltaCCF} --cancer_ctype {params.cancer} "
        "--min_ac_reads {params.min_ac_reads} --clust_dist {params.clust_dist} "

rule CloneGenotype_LongSom:
    input: 
        tsv = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv"
    output:
        tsv = f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
    shell:
        "python {params.scomatic}/scDNAClonesGenotyping/scDNAClonesGenotyping.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} "
        "--nprocs {threads} --min_mq {params.mapq} --tmp_dir {output.tmp}"

rule CloneGenotype_SComatic:
    input: 
        tsv = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam",
        bai = f"{DATA}/bam/scDNA/{{id}}_scDNA.bam.bai",
        barcodes = f"{DATA}/ctypes/scDNA/clones_{{id}}.tsv"
    output:
        tsv = f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
    shell:
        "python {params.scomatic}/scDNAClonesGenotyping/scDNAClonesGenotyping.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} "
        "--nprocs {threads} --min_mq {params.mapq} --tmp_dir {output.tmp}"

rule CloneGenotype_SComaticFromscDNA:
    input: 
        tsv = f"{OUTDIR}/scDNACalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        bam = f"{DATA}/bam/{{id}}.CB.bam",
        bai = f"{DATA}/bam/{{id}}.CB.bam.bai",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv"
    output:
        tsv = f"{OUTDIR}/scDNAValidation/CloneGenotype/SComaticFromscDNA/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/scDNAValidation/CloneGenotype/SComaticFromscDNA/{{id}}/"))
    conda:
        "SComatic"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/CloneGenotype/SComaticFromscDNA",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        ctype=config['CellTypeReannotation']['ctype_column']
    shell:
        "python {params.scomatic}/scDNAClonesGenotyping/scDNAClonesGenotyping.py "
        "--infile {input.tsv} --outfile {params.outdir}/{wildcards.id} "
        "--bam {input.bam} --meta {input.barcodes} --ref {params.hg38} "
        "--nprocs {threads} --min_mq {params.mapq} --tmp_dir {output.tmp} --ctype {params.ctype}"

rule ComparisonSComaticLongSom:
    input:
        scomatic = f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv",
        longsom = f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        scDNACalls = f"{OUTDIR}/scDNACalling/BaseCellCalling/{{id}}.calling.step3.tsv",
        scDNAValidLong = f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv",
        scDNAValidSCom = f"{OUTDIR}/scDNAValidation/CloneGenotype/SComatic/{{id}}.CloneGenotype.tsv",
        scDNA_supp_in_scRNA = f"{OUTDIR}/scDNAValidation/CloneGenotype/SComaticFromscDNA/{{id}}.CloneGenotype.tsv",
    output:
        f"{OUTDIR}/scDNAValidation/Venn3/{{id}}.Venn3.png",
        f"{OUTDIR}/scDNAValidation/Venn3/{{id}}.LongSomScores.tsv",
        f"{OUTDIR}/scDNAValidation/Venn3/{{id}}.SComaticScores.tsv",
    resources:
        time = 120,
        mem_mb = 8000
    params:
        scomatic=SCOMATIC_PATH,
        outdir=f"{OUTDIR}/scDNAValidation/Venn3",
    shell:
        "python {params.scomatic}/ComparisonSComaticLongSom/ComparisonSComaticLongSom.py "
        "--SComatic {input.scomatic} --LongSom {input.longsom} --scDNACalls {input.scDNACalls} "
        "--scDNAValidLong {input.scDNAValidLong} --scDNAValidSCom {input.scDNAValidSCom} "
        " --scDNA_supp_in_scRNA {input.scDNA_supp_in_scRNA} --outfile {params.outdir}/{wildcards.id} "

