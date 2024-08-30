subread = subread-2.0.6-source

include: 'CellTypeReannotation.smk'

INFERCNV_PATH = config['Global']['infercnv']

rule all:
    input:
        expand('featurecount/{id}.counts.formated.txt', id = SAMPLES)

rule split_per_bc:
    input:
        bam = f"{DATA}/bam/{{id}}.bam",
        bai = f"{DATA}/bam/{{id}}.bam.bai",
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv"
    output:
        txt = f"{OUTDIR}/InferCNV/featurecount/input/{{id}}.terminado_splitbc.tx"
    params:
        infercnv = INFERCNV_PATH
    threads:
        32
    resources:
        time = 1200,
        mem_mb=8000
    shell:
        "python {params.infercnv}/split_by_bc.py --cpu {threads} "
        "--id {wildcards.id} --input {input.bam} --barcodes {input.bcs}"

rule featureCount:
	input:
		f"{OUTDIR}/InferCNV/featurecount/input/terminado_splitbc.tx"
	output:
		f"{OUTDIR}/InferCNV/featurecount/output/{{id}}.counts.txt"
    threads:
		32
	resources:
		time = 1200,
        mem_mb=1000
	params:
		anno=config['Global']['isoforms'],	
        subread=config['Global']['subread'],
	shell:
		"{params.subread}/bin/featureCounts -T {threads} -L "
		"-a  {params.anno} -o {output}  barcodes/{wildcards.id}/*.bam "

rule format_featureCount:
	input:
		f"{OUTDIR}/InferCNV/featurecount/output/{{id}}.counts.txt"
	output:
		f"{OUTDIR}/InferCNV/featurecount/output/{{id}}.counts.formated.txt"
	resources:
		time = 1200,
	mem_mb=8000
	run:
		import pandas as pd
		def gene_name_dico(file):
			dico = {}
			names = ['chr','source','type','start','end', '_', 'strand', '__', 'info']
			gencode = pd.read_csv(file, sep = '\t', names = names, skiprows = 5)
			for i,row in gencode.iterrows():
				if row['type'] == 'gene':
					ensg = row['info'].split(';')[0].split('"')[1]
					gene = row['info'].split(';')[2].split('"')[1]
					dico[ensg] = gene
			return dico

		dico = gene_name_dico('/cluster/work/bewi/members/dondia/projects/ovarian_cancer/reference/gencode.v36.annotation.gtf')
		df = pd.read_csv(input[0], comment='#', sep='\t')
		cols = [1,2,3,4,5]
		df.drop(df.columns[cols],axis=1,inplace=True)
		df.columns = [i.split('/')[-1].split('.')[0] for i in df.columns]
		df['Geneid'] = df['Geneid'].map(dico)
		df = df[df['Geneid'].notnull()]
		df = df.groupby('Geneid').sum() #sum redundant genes
		df.to_csv(output[0])

rule ctat_infercnv:
	input:
		counts = f"{OUTDIR}/InferCNV/featurecount/output/{{id}}.counts.formated.txt"
        barcodes = f"{OUTDIR}/CellTypeReannotation/ReannotatedCellTypes/{{id}}.tsv"
	output:
		f"{OUTDIR}/InferCNV/InferCNV/{{id}}/infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.png"
	threads: 10
	resources:
		time = 1200,
		mem_mb=8000
    params:
        order = config['inferCNV']['ordering'],
        infercnv = INFERCNV_PATH
        outdir = f"{OUTDIR}/InferCNV/InferCNV/"
	conda:
		'infercnv'
	shell:
        "Rscript {params.infercnv}/infercnv.R "
        "{input.counts} {input.barcodes} {params.order} "
        "{params.outdir}/{wildcards.id}"
