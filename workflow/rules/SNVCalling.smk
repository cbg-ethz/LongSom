### Rules for calling somatic SNVs 


rule SplitBam:
    input:
        bam=f"{INPUT}/bam/{{id}}.bam",
        bai=f"{INPUT}/bam/{{id}}.bam.bai",
        barcodes="CellTypeReannotation/ReannotatedCellTypes/{id}.tsv" if REANNO else "Barcodes/{id}.tsv",
    output:
        expand("SNVCalling/SplitBam/{{id}}.{celltype}.bam", celltype=['Cancer','Non-Cancer'])
    params:
        script=str(workflow.basedir)+"/scripts/PreProcessing/SplitBamCellTypes.py",
        mapq=config['SNVCalling']['BaseCellCounter']['min_mapping_quality'],
    conda:
        "../envs/SComatic.yaml",
    log:
        "logs/SplitBam/{id}.log",
    benchmark:
        "benchmarks/SplitBam/{id}.benchmark.txt",
    shell:
        r"""
        python {params.script} \
        --bam {input.bam} \
        --meta {input.barcodes} \
        --id {wildcards.id} \
        --outdir SNVCalling/SplitBam \
        --min_MQ {params.mapq}
        """
        
rule BaseCellCounter:
    input:
        bam="SNVCalling/SplitBam/{id}.{celltype}.bam",
        ref=str(workflow.basedir)+config['Reference']['genome'],
    output:
        tsv="SNVCalling/BaseCellCounter/{id}/{id}.{celltype}.tsv",
        tmp=temp(directory("SNVCalling/BaseCellCounter/{id}/temp_{celltype}/"))
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCounter.py",
        chrom=config['SNVCalling']['BaseCellCounter']['chromosomes'],
        mapq=config['SNVCalling']['BaseCellCounter']['min_mapping_quality'],
    conda:
        "../envs/SComatic.yaml"
    threads: 64
    resources:
        mem_mb_per_cpu=1024
    log:
        "logs/BaseCellCounter/{id}.{celltype}.log",
    benchmark:
        "benchmarks/BaseCellCounter/{id}.{celltype}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --bam {input.bam} \
        --ref {input.ref} \
        --chrom {params.chrom} \
        --out_folder SNVCalling/BaseCellCounter/{wildcards.id}/ \
        --nprocs {threads} \
        --min_mq {params.mapq} \
        --tmp_dir {output.tmp}
        """

rule MergeCounts:
    input:
        expand("SNVCalling/BaseCellCounter/{{id}}/{{id}}.{celltype}.tsv", 
            celltype=['Cancer','Non-Cancer'])
    output:
        tsv="SNVCalling/MergeCounts/{id}.BaseCellCounts.AllCellTypes.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/MergeBaseCellCounts.py",
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/MergeCounts/{id}.log",
    benchmark:
        "benchmarks/MergeCounts/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --tsv_folder SNVCalling/BaseCellCounter/{wildcards.id}/ \
        --outfile {output.tsv}
        """

if PON:
    rule BaseCellCalling_step1:
        input: 
            bb="PoN/PoN/BetaBinEstimates.txt",
            tsv="SNVCalling/MergeCounts/{id}.BaseCellCounts.AllCellTypes.tsv",
            ref=str(workflow.basedir)+config['Reference']['genome'],
        output:
            tsv="SNVCalling/BaseCellCalling/{id}.calling.step1.tsv"
        params:
            script=str(workflow.basedir)+"/scripts/SNVCalling/aseCellCalling.step1.py",
            min_cell_types = config['SNVCalling']['BaseCellCalling']['Min_cell_types'],
            min_ac_reads = config['SNVCalling']['BaseCellCalling']['min_ac_reads'],
            min_ac_cells = config['SNVCalling']['BaseCellCalling']['min_ac_cells'],
            alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
            beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
            alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
            beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
        conda:
            "../envs/SComatic.yaml"
        log:
            "logs/BaseCellCalling_step1/{id}.log",
        benchmark:
            "benchmarks/BaseCellCalling_step1/{id}.benchmark.txt"
        shell:
            r"""
            python {params.script} \
            --infile {input.tsv} \
            --ref {input.ref} \
            --outfile SNVCalling/BaseCellCalling/{wildcards.id} \
            --min_cell_types {params.min_cell_types} \
            --min_ac_reads {params.min_ac_reads} \
            --min_ac_cells {params.min_ac_cells} \
            --alpha1 {params.alpha1} \
            --beta1 {params.beta1} \
            --alpha2 {params.alpha2} \
            --beta2 {params.beta2} 
            """

else:
    rule BaseCellCalling_step1:
        input: 
            tsv="SNVCalling/MergeCounts/{id}.BaseCellCounts.AllCellTypes.tsv",
            ref=str(workflow.basedir)+config['Reference']['genome'],
        output:
            tsv="SNVCalling/BaseCellCalling/{id}.calling.step1.tsv"
        params:
            script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCalling.step1.py",
            min_cell_types=config['SNVCalling']['BaseCellCalling']['Min_cell_types'],
            min_ac_reads = config['SNVCalling']['BaseCellCalling']['min_ac_reads'],
            min_ac_cells = config['SNVCalling']['BaseCellCalling']['min_ac_cells'],
            alpha1=config['SNVCalling']['BaseCellCalling']['alpha1'],
            beta1=config['SNVCalling']['BaseCellCalling']['beta1'],
            alpha2=config['SNVCalling']['BaseCellCalling']['alpha2'],
            beta2=config['SNVCalling']['BaseCellCalling']['beta2'],
        conda:
            "../envs/SComatic.yaml"
        log:
            "logs/BaseCellCalling_step1/{id}.log",
        benchmark:
            "benchmarks/BaseCellCalling_step1/{id}.benchmark.txt"
        shell:
            r"""
            python {params.script} \
            --infile {input.tsv} \
            --ref {input.ref} \
            --outfile SNVCalling/BaseCellCalling/{wildcards.id} \
            --min_cell_types {params.min_cell_types} \
            --min_ac_reads {params.min_ac_reads} \
            --min_ac_cells {params.min_ac_cells} \
            --alpha1 {params.alpha1} \
            --beta1 {params.beta1} \
            --alpha2 {params.alpha2} \
            --beta2 {params.beta2}
            """

rule BaseCellCalling_step2:
    input: 
        tsv="SNVCalling/BaseCellCalling/{id}.calling.step1.tsv",
        pon_LR="PoN/PoN/PoN_LR.tsv" if PON else [],
        pon_SR = str(workflow.basedir)+config['Reference']['PoN_SR'],
        RNA_editing=str(workflow.basedir)+config['Reference']['RNA_editing'],
    output:
        tsv="SNVCalling/BaseCellCalling/{id}.calling.step2.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCalling.step2.py",
        min_cell_types = config['SNVCalling']['BaseCellCalling']['Min_cell_types'],
        min_distance = config['SNVCalling']['BaseCellCalling']['min_distance'],
        gnomAD_db=str(workflow.basedir)+config['Reference']['gnomAD_db'],
        max_gnomAD_VAF = config['SNVCalling']['BaseCellCalling']['max_gnomAD_VAF'],
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/BaseCellCalling_step2/{id}.log",
    benchmark:
        "benchmarks/BaseCellCalling_step2/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --infile {input.tsv} \
        --outfile SNVCalling/BaseCellCalling/{wildcards.id} \
        --editing {input.RNA_editing} \
        --pon_SR {input.pon_SR} \
        --pon_LR {input.pon_LR} \
        --gnomAD_db {params.gnomAD_db} \
        --gnomAD_max {params.max_gnomAD_VAF} \
        --min_distance {params.min_distance} 
        """

rule BaseCellCalling_step3:
    input: 
        tsv="SNVCalling/BaseCellCalling/{id}.calling.step2.tsv"
    output:
        tsv="SNVCalling/BaseCellCalling/{id}.calling.step3.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCalling.step3.py",
        deltaVAF=config['SNVCalling']['BaseCellCalling']['deltaVAF'],
        deltaMCF=config['SNVCalling']['BaseCellCalling']['deltaMCF'],
        chrm_conta = config['SNVCalling']['BaseCellCalling']['chrM_contaminant'],
        min_ac_reads = config['SNVCalling']['BaseCellCalling']['min_ac_reads'],
        min_ac_cells = config['SNVCalling']['BaseCellCalling']['min_ac_cells'],
        clust_dist = config['SNVCalling']['BaseCellCalling']['clust_dist'],
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/BaseCellCalling_step3/{id}.log",
    benchmark:
        "benchmarks/BaseCellCalling_step3/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --infile {input.tsv} \
        --outfile SNVCalling/BaseCellCalling/{wildcards.id} \
        --chrM_contaminant {params.chrm_conta} \
        --deltaVAF {params.deltaVAF} \
        --deltaMCF {params.deltaMCF} \
        --min_ac_reads {params.min_ac_reads} \
        --min_ac_cells {params.min_ac_cells} \
        --clust_dist {params.clust_dist}
        """
