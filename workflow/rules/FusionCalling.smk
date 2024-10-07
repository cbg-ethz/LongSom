### Rules to call somatic fusions
 
rule BamToFastq:
    input:
        bam="SNVCalling/SplitBam/{id}.{celltype}.bam"
    output:
        fastq=temp("FusionCalling/{id}/{id}.{celltype}.fastq")
    params:
        script=str(workflow.basedir)+"/scripts/FusionCalling/BamToFastq.py",
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/BamToFastq/{id}.{celltype}.log",
    benchmark:
        "benchmarks/BamToFastq/{id}.{celltype}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --bam {input.bam} \
        --fastq {output.fastq}
        """

rule ConcatFastq:
    input:
        expand("FusionCalling/{{id}}/{{id}}.{celltype}.fastq", celltype=['Cancer','Non-Cancer'])
    output:
        fastq=temp("FusionCalling/{id}/{id}.fastq")
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/ConcatFastq/{id}.log",
    benchmark:
        "benchmarks/ConcatFastq/{id}.benchmark.txt"
    shell:
        r"""
        cat FusionCalling/{wildcards.id}/*.fastq >> {output.fastq}
        """

rule CTATFusion:
    input:
        fastq="FusionCalling/{id}/{id}.fastq"
    output:
        tsv="FusionCalling/{id}/ctat-LR-fusion.fusion_predictions.tsv",
    threads: 16
    resources:
        mem_mb_per_cpu=4096
    container:
        str(workflow.basedir)+"/scripts/FusionCalling/ctat_lr_fusion.v0.13.0.simg"
    log:
        "logs/CTATFusion/{id}.log",
    benchmark:
        "benchmarks/CTATFusion/{id}.benchmark.txt"
    shell:
        r"""
        ctat-LR-fusion \
        -T /output/{input.fastq} \
        --output /output/FusionCalling/{wildcards.id} \
        --genome_lib_dir /ref \
        --prep_reference \
        --CPU {threads} \
        --vis
        """

rule SomaticFusions:
    input:
        fusions="FusionCalling/{id}/ctat-LR-fusion.fusion_predictions.tsv",
        barcodes="CellTypeReannotation/ReannotatedCellTypes/{id}.tsv" if REANNO else "Barcodes/{id}.tsv"
    output:
        fusions="FusionCalling/Somatic/{id}.Fusions.tsv",
        long="FusionCalling/Somatic/{id}.Fusions.SingleCellGenotype.tsv"
    params:
         script=str(workflow.basedir)+"/scripts/FusionCalling/FusionCalling.py",
         deltaMCF=config['FusionCalling']['SomaticFusions']['deltaMCF'],
         min_ac_reads=config['FusionCalling']['SomaticFusions']['min_ac_reads'],
         min_ac_cells=config['FusionCalling']['SomaticFusions']['min_ac_cells'],
         max_MCF_noncancer=config['FusionCalling']['SomaticFusions']['max_MCF_noncancer'],
    conda:
        "../envs/SComatic.yaml"
    log:
        "logs/SomaticFusions/{id}.log",
    benchmark:
        "benchmarks/SomaticFusions/{id}.benchmark.txt"
    shell:
        r"""
        python {params.script} \
        --fusions {input.fusions} \
        --barcodes {input.barcodes} \
        --min_ac_reads {params.min_ac_reads} \
        --min_ac_cells {params.min_ac_cells} \
        --max_MCF_noncancer {params.max_MCF_noncancer} \
        --deltaMCF {params.deltaMCF} \
        --outdir FusionCalling/Somatic/{wildcards.id}
        """