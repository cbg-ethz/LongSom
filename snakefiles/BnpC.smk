import pandas as pd
CTYPES = config['SNVCalling']['celltypes']
OUTDIR=config['Global']['outdir']
OUTBNPC=config['BnpC']['outdir']
DATA=config['Global']['data']
IDS=config['Global']['ids']
BNPC_PATH=config['Global']['bnpc']
SCDNA=config['Run']['scdna']
#include: 'SNVCalling.smk'

rule all_BnpC:
    input:
        expand(f"{OUTDIR}/BnpC/{OUTBNPC}/{{id}}/genoCluster_posterior_mean_raw.pdf",
         id=IDS),
        expand(f"{OUTDIR}/SNVCalling/ClusterMap/{{id}}.ClusterMap.Reannotation.pdf",
         id=IDS),
        expand(f"{OUTDIR}/SNVCalling/Annotations/{{id}}.hg38_multianno.txt", 
         id=IDS),

rule FormatInputBnpC:
    input:
        bin=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.BinaryMatrix.tsv",
        vaf=f"{OUTDIR}/SNVCalling/SingleCellGenotype/{{id}}.VAFMatrix.tsv",
        scDNA=f"{OUTDIR}/scDNAValidation/CloneGenotype/LongSom/{{id}}.CloneGenotype.tsv" if SCDNA else [],
        ctypes=f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv",
    output:
        bin=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.BinaryMatrix.tsv",
        vaf=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.VAFMatrix.tsv",
        scDNA=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.scDNACloneGenotype.tsv",
        ctypes=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.Barcodes.tsv",
    conda:
        "envs/BnpC.yml"
    resources:
        time = 1200,
        mem_mb=2000
    params:
        bnpc = BNPC_PATH,
        outdir=f"{OUTDIR}/BnpC/BnpC_input/",
        min_cells = config['BnpC']['min_cells_per_mut'],
        min_cov = config['BnpC']['min_pos_cov']
    shell:
        "python {params.bnpc}/libs/format_input.py --bin {input.bin} "
        "--vaf {input.vaf} --scDNA {input.scDNA} --ctypes {input.ctypes} "
        "--min_pos_cov {params.min_cov} --min_cells_per_mut {params.min_cells} "
        "--outfile {params.outdir}/{wildcards.id}"


rule BnpC_clones:
    input:
        bin=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.BinaryMatrix.tsv",
        vaf=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.VAFMatrix.tsv",
        scDNA=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.scDNACloneGenotype.tsv",
        ctypes=f"{OUTDIR}/BnpC/BnpC_input/{{id}}.Barcodes.tsv",
    output:
        f"{OUTDIR}/BnpC/{OUTBNPC}/{{id}}/genoCluster_posterior_mean_raw.pdf"
    threads:
        16
    resources:
        time = 1200,
        mem_mb=2000
    params:
        bnpc = BNPC_PATH,
        outdir=f"{OUTDIR}/BnpC/{OUTBNPC}/",
        mcmc_steps = config['BnpC']['mcmc_steps'],
        estimator = config['BnpC']['estimator'],
        dpa = lambda w: config['BnpC']['dpa'][w.get("id")],
        cup = config['BnpC']['cup'],
        eup = config['BnpC']['eup'],
        FP = config['BnpC']['FP'],
        FN = config['BnpC']['FN'],
        pp= config['BnpC']['pp'],
        mut_order= lambda w: config['BnpC']['mut_order'][w.get("id")],
    conda:
        "envs/BnpC.yml"
    shell:
        "python {params.bnpc}/run_BnpC.py {input.bin} -n {threads} "
        "-o {params.outdir}/{wildcards.id} -s {params.mcmc_steps} "
        "-e {params.estimator} -cup {params.cup} -eup {params.eup} "
        "-FP {params.FP} -FN {params.FN} -pp {params.pp} -ap {params.dpa} "
        "--ctypes {input.ctypes} --scDNA {input.scDNA} --mut_order {params.mut_order}"
