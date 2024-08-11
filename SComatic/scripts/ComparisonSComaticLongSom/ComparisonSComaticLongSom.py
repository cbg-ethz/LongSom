import timeit
import argparse
import pandas as pd
import gzip

def clinical_variants(SComatic,LongSom,outfile):
	SComatic = pd.read_csv(SComatic, sep='\t', skiprows=29)
	LongSom = pd.read_csv(LongSom, sep='\t', skiprows=29)

	S = set(SComatic['INDEX'])
	L = set(LongSom['INDEX'])

	shared = [i for i in L if i in S]
	Scomatic_exclusive = [i for i in S if i not in L]
	LongSom_exclusive = [i for i in L if i not in S]

	clinical_variants = pd.read_csv(clinical_anno, sep='\t')
	clinical_variants['INDEX'] = clinical_variants['CHROM'] + ':' + clinical_variants['POS'].astype(str) + ':' + clinical_variants['ALT']
	clinical_variants = clinical_variants[clinical_variants['INDEX'].isin(set(SNVs['INDEX']))]
	clinical_variants.insert(loc=0, column='Patient', value=patient_id)
	clinical_variants = clinical_variants[[i for i in clinical_variants.columns if i not in ['DP','QUAL','MQ']]]
	clinical_variants.to_csv(outfile,sep='\t',index=False)



def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--SComatic', type=str, default=1, help='SComatic base calling file (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--LongSom', type=str, default=1, help='SComatic+CellTypeReannotation base calling file (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--outfile', default = 'Matrix.tsv', help='Out file', required = False)
	parser.add_argument('--id', type=str, help='Patient ID', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	SComatic = args.SComatic
	LongSom = args.LongSom
	out_prefix = args.outfile
	patient_id = args.id

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Create clinical annotation file
	outfile = out_prefix + 'ClinicalAnnotation.tsv'
	clinical_variants(infile,clinical_anno,patient_id,outfile)

	# 2 Create gene annotation file
	outfile = out_prefix + 'GeneAnnotation.tsv'
	gene_annotation(infile,gene_anno,outfile)
	

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



