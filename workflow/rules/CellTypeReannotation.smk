### Rules for cell typeReannotation
### First detects SNVs and fusions
### Then defines the High Confidence Cancer variants (HCCV)
### Then cells areReannotated based on their SNV/fusion HCCV mutation status


### START PREPROCESSING ###
rule RenameCellTypes:
    input:
        barcodes=f"{INPUT}/barcodes/{{id}}.tsv"
    output:
        barcodes="Barcodes/{id}.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/PreProcessing/RenameCellTypes.py",
        cancer_cell_type=config['User']['cancer_cell_type']
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/RenameCellTypes/{id}.log",
    benchmark:
        "benchmarks/RenameCellTypes/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --input {input.barcodes} \
        --output {output.barcodes} \
        --cancer_cell_type {params.cancer_cell_type}
        """

rule SplitBam_Reanno:
    input:
        bam=f"{INPUT}/bam/{{id}}.bam",
        bai=f"{INPUT}/bam/{{id}}.bam.bai",
        barcodes="Barcodes/{id}.tsv"
    output:
        expand("CellTypeReannotation/SplitBam/{{id}}.{celltype}.bam", celltype=['Cancer','Non-Cancer'])
    params:
        script=str(workflow.basedir)+"/scripts/PreProcessing/SplitBamCellTypes.py",
        mapq=config['Reanno']['BaseCellCounter']['min_mapping_quality'],
    conda:
        "../envs/SComatic.yaml",
    log:
        "logs/SplitBam_Reanno/{id}.log",
    benchmark:
        "benchmarks/SplitBam_Reanno/{id}.benchmark.txt",
    shell:
        r"""
        python {params.script} \
        --bam {input.bam} \
        --meta {input.barcodes} \
        --id {wildcards.id} \
        --outdir CellTypeReannotation/SplitBam \
        --min_MQ {params.mapq}
        """
### END PREPROCESSING ###

### START FUSION CALLING ###
rule BamToFastq_Reanno:
    input:
        bam="CellTypeReannotation/SplitBam/{id}.{celltype}.bam"
    output:
        fastq=temp("CellTypeReannotation/FusionCalling/{id}/{id}.{celltype}.fastq")
    params:
        script=str(workflow.basedir)+"/scripts/FusionCalling/BamToFastq.py",
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/BamToFastq_Reanno/{id}.{celltype}.log",
    benchmark:
        "benchmarks/BamToFastq_Reanno/{id}.{celltype}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --bam {input.bam} \
        --fastq {output.fastq}
        """

rule ConcatFastq_Reanno:
    input:
        expand("CellTypeReannotation/FusionCalling/{{id}}/{{id}}.{celltype}.fastq", celltype=['Cancer','Non-Cancer'])
    output:
        fastq=temp("CellTypeReannotation/FusionCalling/{id}/{id}.fastq")
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/ConcatFastq_Reanno/{id}.log",
    benchmark:
        "benchmarks/ConcatFastq_Reanno/{id}.benchmark.txt"
    shell:
        r"""
        cat CellTypeReannotation/FusionCalling/{wildcards.id}/*.fastq >> {output.fastq}
        """

rule CTATFusion_Reanno:
    input:
        fastq="CellTypeReannotation/FusionCalling/{id}/{id}.fastq",
    output:
        tsv="CellTypeReannotation/FusionCalling/{id}/ctat-LR-fusion.fusion_predictions.tsv",
    threads: 16
    resources:
        mem_mb_per_cpu=4096
    container:
        str(workflow.basedir)+"/scripts/FusionCalling/ctat_lr_fusion.v0.13.0.simg"
    log:
        "logs/CTATFusion_Reanno/{id}.log",
    benchmark:
        "benchmarks/CTATFusion_Reanno/{id}.benchmark.txt"
    shell:
        r"""
        ctat-LR-fusion \
        -T /output/{input.fastq} \
        --output /output/CellTypeReannotation/FusionCalling/{wildcards.id} \
        --genome_lib_dir /ref \
        --prep_reference \
        --CPU {threads} \
        --vis
        """
### END FUSION CALLING ###

### START SNV CALLING ###
rule BaseCellCounter_Reanno:
    input:
        bam="CellTypeReannotation/SplitBam/{id}.{celltype}.bam",
        ref=str(workflow.basedir)+config['Reference']['genome'],
    output:
        tsv="CellTypeReannotation/BaseCellCounter/{id}/{id}.{celltype}.tsv",
        tmp=temp(directory("CellTypeReannotation/BaseCellCounter/{id}/temp_{celltype}/"))
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCounter.py",
        chrom=config['Reanno']['BaseCellCounter']['chromosomes'],
        mapq=config['Reanno']['BaseCellCounter']['min_mapping_quality'],
    conda:
        "../envs/SComatic.yaml"
    threads: 64
    resources:
        mem_mb_per_cpu=1024
    log:
        "logs/BaseCellCounter_Reanno/{id}.{celltype}.log",
    benchmark:
        "benchmarks/BaseCellCounter_Reanno/{id}.{celltype}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --bam {input.bam} \
        --ref {input.ref} \
        --chrom {params.chrom} \
        --out_folder CellTypeReannotation/BaseCellCounter/{wildcards.id}/ \
        --nprocs {threads} \
        --min_mq {params.mapq} \
        --tmp_dir {output.tmp}
        """

rule MergeCounts_Reanno:
    input:
        expand("CellTypeReannotation/BaseCellCounter/{{id}}/{{id}}.{celltype}.tsv", 
            celltype=['Cancer','Non-Cancer'])
    output:
        tsv="CellTypeReannotation/MergeCounts/{id}.BaseCellCounts.AllCellTypes.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/MergeBaseCellCounts.py",
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/MergeCounts_Reanno/{id}.log",
    benchmark:
        "benchmarks/MergeCounts_Reanno/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --tsv_folder CellTypeReannotation/BaseCellCounter/{wildcards.id}/ \
        --outfile {output.tsv}
        """
if PON:
	rule BaseCellCalling_step1_Reanno:
		input: 
			bb="PoN/PoN/BetaBinEstimates.txt",
			tsv="CellTypeReannotation/MergeCounts/{id}.BaseCellCounts.AllCellTypes.tsv",
			ref= str(workflow.basedir)+config['Reference']['genome'],
		output:
			tsv="CellTypeReannotation/BaseCellCalling/{id}.calling.step1.tsv"
		params:
			script=str(workflow.basedir)+"/scripts/SNVCalling/aseCellCalling.step1.py",
			min_cell_types = config['Reanno']['BaseCellCalling']['Min_cell_types'],
			alpha1 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha1'),
			beta1 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta1'),
			alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
			beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
		conda:
			"../envs/SComatic.yaml"
		log:
			"logs/BaseCellCalling_step1_Reanno/{id}.log",
		benchmark:
			"benchmarks/BaseCellCalling_step1_Reanno/{id}.benchmark.txt"
		shell:
			r"""
			python {params.script} \
			--infile {input.tsv} \
			--ref {input.ref} \
			--outfile CellTypeReannotation/BaseCellCalling/{wildcards.id} \
			--min_cell_types {params.min_cell_types} \
			--alpha1 {params.alpha1} \
			--beta1 {params.beta1} \
			--alpha2 {params.alpha2} \
			--beta2 {params.beta2} 
			"""

else:
	rule BaseCellCalling_step1_Reanno:
		input: 
			tsv="CellTypeReannotation/MergeCounts/{id}.BaseCellCounts.AllCellTypes.tsv",
			ref=str(workflow.basedir)+config['Reference']['genome']
		output:
			tsv="CellTypeReannotation/BaseCellCalling/{id}.calling.step1.tsv"
		params:
			script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCalling.step1.py",
			min_cell_types=config['Reanno']['BaseCellCalling']['Min_cell_types'],
			alpha1=config['Reanno']['BaseCellCalling']['alpha1'],
			beta1=config['Reanno']['BaseCellCalling']['beta1'],
			alpha2=config['Reanno']['BaseCellCalling']['alpha2'],
			beta2=config['Reanno']['BaseCellCalling']['beta2'],
		conda:
			"../envs/SComatic.yaml"
		log:
			"logs/BaseCellCalling_step1_Reanno/{id}.log",
		benchmark:
			"benchmarks/BaseCellCalling_step1_Reanno/{id}.benchmark.txt"
		shell:
			r"""
			python {params.script} \
			--infile {input.tsv} \
			--ref {input.ref} \
			--outfile CellTypeReannotation/BaseCellCalling/{wildcards.id} \
			--min_cell_types {params.min_cell_types} \
			--alpha1 {params.alpha1} \
			--beta1 {params.beta1} \
			--alpha2 {params.alpha2} \
			--beta2 {params.beta2}
			"""

rule BaseCellCalling_step2_Reanno:
    input: 
        tsv="CellTypeReannotation/BaseCellCalling/{id}.calling.step1.tsv",
        pon_LR="PoN/PoN/PoN_LR.tsv" if PON else [],
    output:
        tsv="CellTypeReannotation/BaseCellCalling/{id}.calling.step2.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/SNVCalling/BaseCellCalling.step2.py",
        min_cell_types=config['Reanno']['BaseCellCalling']['Min_cell_types'],
        min_distance=config['Reanno']['BaseCellCalling']['min_distance'],
        gnomAD_db=str(workflow.basedir)+config['Reference']['gnomAD_db'],
        pon_SR = str(workflow.basedir)+config['Reference']['PoN_SR'],
        RNA_editing = str(workflow.basedir)+config['Reference']['RNA_editing'],
        max_gnomAD_VAF = config['Reanno']['BaseCellCalling']['max_gnomAD_VAF'],
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/BaseCellCalling_step2_Reanno/{id}.log",
    benchmark:
        "benchmarks/BaseCellCalling_step2_Reanno/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --infile {input.tsv} \
        --outfile CellTypeReannotation/BaseCellCalling/{wildcards.id} \
        --editing {params.RNA_editing} \
        --pon_SR {params.pon_SR} \
        --pon_LR {input.pon_LR} \
        --gnomAD_db {params.gnomAD_db} \
        --gnomAD_max {params.max_gnomAD_VAF} \
        --min_distance {params.min_distance} 
        """
### END SNV CALLING ###

### START HCCV CALLING ###
rule HighConfidenceCancerVariants:
    input: 
        tsv="CellTypeReannotation/BaseCellCalling/{id}.calling.step2.tsv"
    output:
        tsv="CellTypeReannotation/HCCV/{id}.HCCV.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/CellTypeReannotation/HighConfidenceCancerVariants.py",
        min_dp=config['Reanno']['HCCV']['min_depth'],
        deltaVAF=config['Reanno']['HCCV']['deltaVAF'],
        deltaMCF=config['Reanno']['HCCV']['deltaMCF'],
        clust_dist = config['Reanno']['HCCV']['clust_dist'],
    conda:
        "../envs/SComatic.yaml"
    threads: 8
    resources:
        mem_mb_per_cpu=1024
    log:
        "logs/HighConfidenceCancerVariants/{id}.log",
    benchmark:
        "benchmarks/HighConfidenceCancerVariants/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --SNVs {input.tsv} \
        --outfile CellTypeReannotation/HCCV/{wildcards.id} \
        --min_dp {params.min_dp} \
        --deltaVAF {params.deltaVAF} \
        --deltaMCF {params.deltaMCF} \
        --clust_dist {params.clust_dist}
        """

if PON:
    rule HCCVSingleCellGenotype:
        input: 
            tsv="CellTypeReannotation/HCCV/{id}.HCCV.tsv",
            bam=f"{INPUT}/bam/{{id}}.bam",
            barcodes="Barcodes/{id}.tsv",
            bb="PoN/PoN/BetaBinEstimates.txt",
            ref=str(workflow.basedir)+config['Reference']['genome'],
        output:
            tsv="CellTypeReannotation/HCCV/{id}.SNVs.SingleCellGenotype.tsv",
            tmp=temp(directory("CellTypeReannotation/HCCV/{id}/"))
        params:
            script=str(workflow.basedir)+"/scripts/CellTypeReannotation/HCCVSingleCellGenotype.py",
            alt_flag= config['Reanno']['HCCV']['alt_flag'],
            mapq=config['Reanno']['BaseCellCounter']['min_mapping_quality'],
            alpha2 = lambda w, input: get_BetaBinEstimates(input.bb, 'alpha2'),
            beta2 = lambda w, input: get_BetaBinEstimates(input.bb, 'beta2'),
            pval = config['Reanno']['HCCV']['pvalue'],
            chrm_conta = config['Reanno']['HCCV']['chrM_contaminant']
        conda:
            "../envs/SComatic.yaml"
        threads: 32
        resources:
            mem_mb_per_cpu=1024
        log:
            "logs/HCCVSingleCellGenotype/{id}.log",
        benchmark:
            "benchmarks/HCCVSingleCellGenotype/{id}.benchmark.txt"
        shell:
            r"""
            python {params.script} \
            --bam {input.bam} \
            --infile {input.tsv} \
            --ref {input.ref} \
            --outfile {output.tsv} \
            --meta {input.barcodes} \
            --alt_flag {params.alt_flag} \
            --nprocs {threads} \
            --min_mq {params.mapq} \
            --pvalue {params.pval} \
            --alpha2 {params.alpha2} \
            --beta2 {params.beta2} \
            --chrM_contaminant {params.chrm_conta} \
            --tmp_dir {output.tmp}
            """
else:
    rule HCCVSingleCellGenotype:
        input: 
            tsv="CellTypeReannotation/HCCV/{id}.HCCV.tsv",
            bam=f"{INPUT}/bam/{{id}}.bam",
            barcodes="Barcodes/{id}.tsv",
            ref=str(workflow.basedir)+config['Reference']['genome']
        output:
            tsv="CellTypeReannotation/HCCV/{id}.SNVs.SingleCellGenotype.tsv",
            tmp=temp(directory("CellTypeReannotation/HCCV/{id}/"))
        params:
            script=str(workflow.basedir)+"/scripts/CellTypeReannotation/HCCVSingleCellGenotype.py",
            alt_flag= config['Reanno']['HCCV']['alt_flag'],
            mapq=config['Reanno']['BaseCellCounter']['min_mapping_quality'],
            alpha2 = config['Reanno']['BaseCellCalling']['alpha2'],
            beta2 = config['Reanno']['BaseCellCalling']['beta2'],
            pval = config['Reanno']['HCCV']['pvalue'],
            chrm_conta = config['Reanno']['HCCV']['chrM_contaminant']
        conda:
            "../envs/SComatic.yaml"
        threads: 32
        resources:
            mem_mb_per_cpu=1024
        log:
            "logs/HCCVSingleCellGenotype/{id}.log",
        benchmark:
            "benchmarks/HCCVSingleCellGenotype/{id}.benchmark.txt"
        shell:
            r"""
            python {params.script} \
            --bam {input.bam} \
            --infile {input.tsv} \
            --ref {input.ref} \
            --outfile {output.tsv} \
            --meta {input.barcodes} \
            --alt_flag {params.alt_flag} \
            --nprocs {threads} \
            --min_mq {params.mapq} \
            --pvalue {params.pval} \
            --alpha2 {params.alpha2} \
            --beta2 {params.beta2} \
            --chrM_contaminant {params.chrm_conta} \
            --tmp_dir {output.tmp}
            """

rule HCCVFusions:
    input:
        fusions="CellTypeReannotation/FusionCalling/{id}/ctat-LR-fusion.fusion_predictions.tsv",
        barcodes="Barcodes/{id}.tsv"
    output:
        fusions="CellTypeReannotation/HCCV/{id}.Fusions.tsv",
        long="CellTypeReannotation/HCCV/{id}.Fusions.SingleCellGenotype.tsv"
    params:
         script=str(workflow.basedir)+"/scripts/FusionCalling/FusionCalling.py",
         deltaMCF=config['FusionCalling']['SomaticFusions']['deltaMCF'],
         min_ac_reads=config['FusionCalling']['SomaticFusions']['min_ac_reads'],
         min_ac_cells=config['FusionCalling']['SomaticFusions']['min_ac_cells'],
         max_MCF_noncancer=config['FusionCalling']['SomaticFusions']['max_MCF_noncancer'],
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/HCCV_Fusion/{id}.log",
    benchmark:
        "benchmarks/HCCV_Fusion/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --fusions {input.fusions} \
        --barcodes {input.barcodes} \
        --min_ac_reads {params.min_ac_reads} \
        --min_ac_cells {params.min_ac_cells} \
        --max_MCF_noncancer {params.max_MCF_noncancer} \
        --deltaMCF {params.deltaMCF} \
        --outdir CellTypeReannotation/HCCV/{wildcards.id} 
        """
### END HCCV CALLING ###

### REANNOTATION ###
rule CellTypeReannotation:
    input:
        SNVs="CellTypeReannotation/HCCV/{id}.SNVs.SingleCellGenotype.tsv",
        fusions="CellTypeReannotation/HCCV/{id}.Fusions.SingleCellGenotype.tsv",
        barcodes="Barcodes/{id}.tsv"
    output:
        barcodes="CellTypeReannotation/ReannotatedCellTypes/{id}.tsv"
    params:
        script=str(workflow.basedir)+"/scripts/CellTypeReannotation/CellTypeReannotation.py",
        min_variants = config['Reanno']['Reannotation']['min_variants'],
        min_frac = config['Reanno']['Reannotation']['min_fraction']
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/CellTypeReannotation/{id}.log",
    benchmark:
        "benchmarks/CellTypeReannotation/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --SNVs {input.SNVs} \
        --fusions {input.fusions} \
        --outfile {output.barcodes} \
        --meta {input.barcodes} \
        --min_variants {params.min_variants} \
        --min_frac {params.min_frac}
        """

