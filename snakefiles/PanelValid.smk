OUTDIR=config['Global']['outdir'] 
IDS=config['Global']['ids']
SCOMATIC_PATH=config['Global']['scomatic']
QC_PATH=config['Global']['qc']
DATA=config['Global']['data']

import pandas as pd
def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]

rule all:
    input:
        expand(f"{OUTDIR}/PanelValidation/{{id}}.CloneGenotype.tsv", id=IDS),
        expand(f"{OUTDIR}/PanelValidation/PanelValidStats.tsv", id=IDS),
        expand(f"{OUTDIR}/SNVCalling/Annotations/{{id}}.hg38_multianno.txt", id=IDS),

rule Panel_supp_in_scRNA:
    input: 
        tsv = f"{OUTDIR}/PanelValidation/{{id}}.Panel.tsv",
        bam = f"{DATA}/bam/{{id}}.bam",
        bai = f"{DATA}/bam/{{id}}.bam.bai",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv",
        bb = f"{OUTDIR}/PoN/PoN/BetaBinEstimates.txt",
    output:
        tsv = f"{OUTDIR}/PanelValidation/{{id}}.CloneGenotype.tsv",
        tmp=temp(directory(f"{OUTDIR}/PanelValidation/{{id}}/"))
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 8000
    params:
        qc=QC_PATH,
        outdir=f"{OUTDIR}/PanelValidation",
        hg38=config['Global']['genome'],
        alt_flag= config['SComatic']['SingleCellGenotype']['alt_flag'],
        mapq=config['scDNA']['BaseCellCounter']['min_mapping_quality'],
        ctype="Cell_type",
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
        "--nprocs 1 --tmp_dir {output.tmp} --ctype {params.ctype}"

rule PanelValid:
    input: 
        expand(f"{OUTDIR}/PanelValidation/{{id}}.Panel.tsv", id=IDS),
        expand(f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/{{id}}.calling.step3.tsv", id=IDS),
        expand(f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv", id=IDS),
    output:
        f"{OUTDIR}/PanelValidation/PanelValidStats.tsv"
    conda:
        "envs/SComatic.yml"
    threads:
        32
    resources:
        time = 120,
        mem_mb = 1000
    params:
        qc=QC_PATH,
        panel=f"{OUTDIR}/PanelValidation",
        longsom=f"{OUTDIR}/SNVCalling/BaseCellCalling/",
        scomatic=f"{OUTDIR}/CellTypeReannotation/BaseCellCalling/",
    shell:
        "python {params.qc}/PanelValidation/PanelValidation.py "
        "--panel_dir {params.panel} --longsom_dir {params.longsom} "
        "--scomatic_dir {params.scomatic} --outfile {output} "

rule input_annovar:
    input:
        f"{OUTDIR}/SNVCalling/BaseCellCalling/{{id}}.calling.step3.tsv",
    output:
        f"{OUTDIR}/SNVCalling/Annotations/{{id}}.avinput" 
    shell:
        """grep -v '#' {input} |  tr '\t' '-' | """
        """awk -F'-' -v OFS='\t' '{{print $1,$2,$3,$4,$5,$0}}' | """
        """sed 's/A,A/A/g' | sed 's/C,C/C/g' | sed 's/G,G/G/g' | sed 's/T,T/T/g' """
        """> {output} """

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
        "-protocol refGene,exac03,avsnp147,dbnsfp47a,clinvar_20240611,dbscsnv11,cosmic70 "
        "-operation gx,f,f,f,f,f,f "

# rule clinical_info:
#     input:
#         expand(f"{OUTDIR}/SNVCalling/Annotations/{{id}}.hg38_multianno.txt", id=IDS)
#     output:
#         f"{OUTDIR}/SNVCalling/Annotations/Clinical_Variants.tsv