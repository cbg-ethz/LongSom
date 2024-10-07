### Rules for establishing the genotype of each cells (i.e. SNVs and fusions per cell)
### Then building a cell-variant matrix and clustering cells based on variants (using BnpC)
if PON:
    rule SingleCellGenotype:
        input: 
            tsv="SNVCalling/BaseCellCalling/{id}.calling.step3.tsv",
            bam=f"{INPUT}/bam/{{id}}.bam",
            barcodes="CellTypeReannotation/ReannotatedCellTypes/{id}.tsv",
            bb="PoN/PoN/BetaBinEstimates.txt",
            fusions="FusionCalling/Somatic/{id}.Fusions.SingleCellGenotype.tsv" if CTATFUSION else [],
            ref=str(workflow.basedir)+config['Reference']['genome'],
        output:
            tsv="CellClustering/SingleCellGenotype/{id}.SingleCellGenotype.tsv",
            dp="CellClustering/SingleCellGenotype/{id}.DpMatrix.tsv",
            alt="CellClustering/SingleCellGenotype/{id}.AltMatrix.tsv",
            vaf="CellClustering/SingleCellGenotype/{id}.VAFMatrix.tsv",
            bin="CellClustering/SingleCellGenotype/{id}.BinaryMatrix.tsv",
            tmp=temp(directory("CellClustering/SingleCellGenotype/{id}/"))
        params:
            script=str(workflow.basedir)+"/scripts/CellClustering/SingleCellGenotype.py",
            alt_flag= config['CellClust']['SingleCellGenotype']['alt_flag'],
            mapq=config['SNVCalling']['BaseCellCounter']['min_mapping_quality'],
            alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
            beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
            pval = config['CellClust']['SingleCellGenotype']['pvalue'],
            chrm_conta = config['SNVCallling']['BaseCellCalling']['chrM_contaminant'],
        conda:
            "../envs/SComatic.yaml"
        threads: 32
        resources:
            mem_mb_per_cpu=1024
        log:
            "logs/SingleCellGenotype/{id}.log",
        benchmark:
            "benchmarks/SingleCellGenotype/{id}.benchmark.txt"
        shell:
            r"""
            python {params.script} \
            --infile {input.tsv} \
            --outfile CellClustering/SingleCellGenotype/{wildcards.id} \
            --bam {input.bam} \
            --meta {input.barcodes} \
            --ref {input.ref} \
            --fusions {input.fusions} \
            --nprocs {threads} \
            --min_mq {params.mapq} \
            --pvalue {params.pval} \
            --alpha2 {params.alpha2} \
            --beta2 {params.beta2} \
            --alt_flag {params.alt_flag} \
            --chrM_contaminant {params.chrm_conta} \
            --tmp_dir {output.tmp} 
            """
else:
    rule SingleCellGenotype:
        input: 
            tsv="SNVCalling/BaseCellCalling/{id}.calling.step3.tsv",
            bam=f"{INPUT}/bam/{{id}}.bam",
            barcodes="CellTypeReannotation/ReannotatedCellTypes/{id}.tsv",
            fusions="FusionCalling/Somatic/{id}.Fusions.tsv" if CTATFUSION else [],
            ref=str(workflow.basedir)+config['Reference']['genome'],
        output:
            tsv="CellClustering/SingleCellGenotype/{id}.SingleCellGenotype.tsv",
            dp="CellClustering/SingleCellGenotype/{id}.DpMatrix.tsv",
            alt="CellClustering/SingleCellGenotype/{id}.AltMatrix.tsv",
            vaf="CellClustering/SingleCellGenotype/{id}.VAFMatrix.tsv",
            bin="CellClustering/SingleCellGenotype/{id}.BinaryMatrix.tsv",
            tmp=temp(directory("CellClustering/SingleCellGenotype/{id}/"))
        params:
            script=str(workflow.basedir)+"/scripts/CellClustering/SingleCellGenotype.py",
            alt_flag= config['CellClust']['SingleCellGenotype']['alt_flag'],
            mapq=config['SNVCalling']['BaseCellCounter']['min_mapping_quality'],
            alpha2=config['SNVCalling']['BaseCellCalling']['alpha2'],
            beta2=config['SNVCalling']['BaseCellCalling']['beta2'],
            pval=config['CellClust']['SingleCellGenotype']['pvalue'],
            chrm_conta=config['SNVCalling']['BaseCellCalling']['chrM_contaminant'],
        conda:
            "../envs/SComatic.yaml"
        threads: 32
        resources:
            mem_mb_per_cpu=1024
        log:
            "logs/SingleCellGenotype/{id}.log",
        benchmark:
            "benchmarks/SingleCellGenotype/{id}.benchmark.txt"
        shell:
            r"""
            python {params.script} \
            --infile {input.tsv} \
            --outfile CellClustering/SingleCellGenotype/{wildcards.id} \
            --bam {input.bam} \
            --meta {input.barcodes} \
            --ref {input.ref} \
            --fusions {input.fusions} \
            --nprocs {threads} \
            --min_mq {params.mapq} \
            --pvalue {params.pval} \
            --alpha2 {params.alpha2} \
            --beta2 {params.beta2} \
            --alt_flag {params.alt_flag} \
            --chrM_contaminant {params.chrm_conta} \
            --tmp_dir {output.tmp} 
            """

rule FormatInputBnpC:
    input:
        bin="CellClustering/SingleCellGenotype/{id}.BinaryMatrix.tsv",
        vaf="CellClustering/SingleCellGenotype/{id}.VAFMatrix.tsv",
        ctypes="CellTypeReannotation/ReannotatedCellTypes/{id}.tsv",
    output:
        bin="CellClustering/BnpC_input/{id}.BinaryMatrix.tsv",
        vaf="CellClustering/BnpC_input/{id}.VAFMatrix.tsv",
        ctypes="CellClustering/BnpC_input/{id}.Barcodes.tsv",
    params:
        script=str(workflow.basedir)+"/scripts/CellClustering/FormatInputBnpC.py",
        min_cells=config['CellClust']['FormatInput']['min_cells_per_mut'],
        min_cov=config['CellClust']['FormatInput']['min_pos_cov']
    conda:
        "../envs/BnpC.yaml"
    log:
        "logs/FormatInputBnpC/{id}.log",
    benchmark:
        "benchmarks/FormatInputBnpC/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --bin {input.bin} \
        --vaf {input.vaf} \
        --ctypes {input.ctypes} \
        --min_pos_cov {params.min_cov} \
        --min_cells_per_mut {params.min_cells} \
        --outfile CellClustering/BnpC_input//{wildcards.id} 
        """

rule BnpC_clustering:
    input:
        bin="CellClustering/BnpC_input/{id}.BinaryMatrix.tsv",
        vaf="CellClustering/BnpC_input/{id}.VAFMatrix.tsv",
        ctypes="CellClustering/BnpC_input/{id}.Barcodes.tsv",
    output:
        pdf="CellClustering/BnpC_output/{id}/genoCluster_posterior_mean_raw.pdf"
    params:
        script=str(workflow.basedir)+"/scripts/CellClustering/run_BnpC.py",
        mcmc_steps = config['CellClust']['BnpC']['mcmc_steps'],
        estimator = config['CellClust']['BnpC']['estimator'],
        dpa = config['CellClust']['BnpC']['dpa'],
        cup = config['CellClust']['BnpC']['cup'],
        eup = config['CellClust']['BnpC']['eup'],
        FP = config['CellClust']['BnpC']['FP'],
        FN = config['CellClust']['BnpC']['FN'],
        pp= config['CellClust']['BnpC']['pp'],
    conda:
        "../envs/BnpC.yaml"
    threads: 16
    resources:
        mem_mb_per_cpu=1024
    log:
        "logs/BnpC/{id}.log",
    benchmark:
        "benchmarks/BnpC/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        {input.bin} \
        -n {threads} \
        -o CellClustering/BnpC_output/{wildcards.id} \
        -s {params.mcmc_steps} \
        -e {params.estimator} \
        -cup {params.cup} \
        -eup {params.eup} \
        -FP {params.FP} \
        -FN {params.FN} \
        -pp {params.pp} \
        -ap {params.dpa} \
        --ctypes {input.ctypes} 
        """
