import timeit
import argparse
import pandas as pd
from itertools import chain
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles

def Venn2(scDNAValid_L,scDNAValid_S,out_prefix):
	# Defining mutation sets: without fusions, and only positions with suficient cov in RNA (scDNA cov is ensured by default)
	A = [i for i in scDNAValid_S['INDEX'] if '--' not in i]
	B = [i for i in scDNAValid_L['INDEX'] if '--' not in i]
	Ab = len([i for i in A if i not in B])
	aB = len([i for i in B if i not in A])
	AB = len([i for i in A if i in B])
	plt.figure(figsize=(4,4))
	v = venn2(subsets=(Ab,aB,AB),set_labels = ('SComatic','LongSom'))

	plt.savefig(out_prefix + '.Venn2.png', dpi=600)
	plt.close()
	

def Venn3(dp,scDNAValid_L,scDNAValid_S,scDNACalls,out_prefix):

	d = {'LongSom':{},'SComatic': {}}
	d = genotyping(scDNAValid_L,dp,d,'LongSom')
	d = genotyping(scDNAValid_S,dp,d,'SComatic')

	# Defining mutation sets: without fusions, and only positions with suficient cov in RNA (scDNA cov is ensured by default)
	A = [i for v in d['SComatic'].values() for i in v ]
	B = [i for v in d['LongSom'].values() for i in v ]
	C = list(scDNACalls['INDEX']) + d['LongSom']['PASS'] + d['LongSom']['SUPP'] + d['SComatic']['PASS'] + d['SComatic']['SUPP']

	Abc = len([i for i in set(A) if i not in set(B+C)])
	aBc = len([i for i in set(B) if i not in set(A+C)])
	abC = len([i for i in set(C) if i not in set(A+B)])
	ABc = len([i for i in set(A) if i in set(B) and i not in set(C)])
	aBC = len([i for i in set(C) if i in set(B) and i not in set(A)])
	AbC = len([i for i in set(A) if i in set(C) and i not in set(B)])
	ABC = len([i for i in set(A) if i in set(C) and i in set(B)])

	plt.figure(figsize=(4,4))
	v = venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),set_labels = ('SComatic','LongSom','scWGS'))

	plt.savefig(out_prefix + '.Venn3.png', dpi=600)

	return A,B,C,d

def scDNASupport(A,B,d,dp,scDNAValid_L,scDNAValid_S,sampleid,out_prefix):
	Ltot = list(scDNAValid_L['INDEX'])
	Stot = list(scDNAValid_S['INDEX'])
	print('Ltot', len(Ltot))
	print('Stot', len(Stot))
	L = [i for i in Ltot if i not in Stot]
	S = [i for i in Stot if i not in Ltot]
	C = [i for i in Ltot if i in Stot]
	print('L', len(L))
	print('S', len(S))
	print('C', len(C))
	d = {'LongSom':{},'SComatic': {}, 'Common':{}}
	d = genotyping(scDNAValid_L[scDNAValid_L['INDEX'].isin(L)],dp,d,'LongSom')
	d = genotyping(scDNAValid_S[scDNAValid_S['INDEX'].isin(S)],dp,d,'SComatic')
	d = genotyping(scDNAValid_L[scDNAValid_L['INDEX'].isin(C)],dp,d,'Common')

	with open(out_prefix + '.scDNASupport.tsv','w') as f:
		TITLE = ['SampleID','Method','#Muts','PASS','GERM','SUPP','NOTUM_HIGHCOV','FRAC_SUP','FRAC_GERM']
		f.write('\t'.join(TITLE)+'\n')
		ComputeSupport(f,sampleid,d,'LongSom')
		ComputeSupport(f,sampleid,d,'SComatic')
		ComputeSupport(f,sampleid,d,'Common')


def F1Score(A,B,C,d,RNASupp,sampleid,out_prefix):
	# Defining Validation Set: SNVs detected in scDNA with at least 1 mutated read in scRNA 
	scDNAValidSet = list(RNASupp[(RNASupp['HGSOC_VAF']>0) | (RNASupp['Non-HGSOC_VAF']>0)]['INDEX'])
	
	with open(out_prefix + '.ValidationSet.tsv','w') as f:
		f.write('ValidationSet\n')
		for i in scDNAValidSet:
			f.write(i + '\n')

	with open(out_prefix + '.F1Scores.tsv','w') as f:
		f.write('SampleID\tMethod\tPrecision\tSensitivity\tF1\n')

		### LongSom Precision:
		# True positives re. precision are SNVs detected by LongSom and detected in scDNA 
		# or with at least 1 supporting read in scDNA and no reads supporting any other alternative allele
		TP_prec = len(
				set( [i for i in set(B) if i in set(scDNAValidSet)] + d['LongSom']['PASS'] + d['LongSom']['SUPP'] )
				) 

		# False positives are SNVs detected in scDNA but not in LongSom, and SNVs detected in LongSom but germlines in scDNA
		FP = len( 
				set( [i for i in set(B) if i not in set(scDNAValidSet)] + d['LongSom']['GERM'] ) 
			)
		Precision = TP_prec / (TP_prec + FP)

		# LongSom Sensitivity:
		# True positives re. sensitivity are SNVs detected by LongSom 
		TP_sens = len([i for i in set(B) if i in set(scDNAValidSet)])
		# False negatives are detected in scDNA, but not by LongSom
		FN = len([i for i in set(scDNAValidSet) if i not in set(B)])
		Sensitivity = TP_sens / (TP_sens + FN)

		# LongSom F1 score:
		F1 = (2*TP_prec) / ((2*TP_prec) + FP + FN) 
		
		f.write(f'{sampleid}\tLongSom\t{Precision}\t{Sensitivity}\t{F1}\n')

		### SComatic Precision:
		TP_prec = len(
				set( [i for i in set(A) if i in set(scDNAValidSet)] + d['SComatic']['PASS'] + d['SComatic']['SUPP'] )
				) 
		FP = len( 
				set( [i for i in set(A) if i not in set(scDNAValidSet)] + d['SComatic']['GERM'] ) 
			)
		Precision = TP_prec / (TP_prec + FP)

		# LongSom Sensitivity:
		TP_sens = len([i for i in set(A) if i in set(scDNAValidSet)])
		FN = len([i for i in set(scDNAValidSet) if i not in set(A)])
		Sensitivity = TP_sens / (TP_sens + FN)

		# SComatic F1 score:
		F1 = (2*TP_prec) / ((2*TP_prec) + FP + FN) 

		f.write(f'{sampleid}\tSComatic\t{Precision}\t{Sensitivity}\t{F1}\n')

def genotyping(geno, dp, d, key):
	#Filtering fusions
	geno = geno[~geno['INDEX'].str.contains('--')]
	PASS = list(geno[(geno['Clone_Tum_MutatedStatus']=='PASS') & (geno['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX'])
	GERM = [i for i in geno[geno['Clone_NonTum_VAF']>0]['INDEX'] if i not in PASS]
	SUPP = [i for i in geno[geno['Clone_Tum_VAF']>0]['INDEX'] if i not in PASS+GERM]
	NOTUM_HIGHCOV = [i for i in geno[geno['Clone_Tum_DP']>=dp]['INDEX'] if i not in  PASS+GERM+SUPP]

	d[key]['PASS'] = PASS
	d[key]['GERM'] = GERM
	d[key]['SUPP'] = SUPP
	d[key]['NOTUM_HIGHCOV'] = NOTUM_HIGHCOV

	return d

def ComputeSupport(f,sampleid,d,key):

	PASS = len(d[key]['PASS'])
	GERM = len(d[key]['GERM'])
	SUPP = len(d[key]['SUPP'])
	NOTUM_HIGHCOV = len(d[key]['NOTUM_HIGHCOV'])

	if sum([PASS,GERM,SUPP,NOTUM_HIGHCOV]) > 0:
		FRAC_SUP = (PASS+SUPP)/sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
		FRAC_GERM = (GERM)/sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
	else:
		FRAC_SUP = 0
		FRAC_GERM = 0
	LEN = sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
	LIST = [sampleid,key,LEN,PASS,GERM,SUPP,NOTUM_HIGHCOV,round(FRAC_SUP,2),round(FRAC_GERM,2)]
	f.write('\t'.join([str(i) for i in LIST])+'\n')

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--dp', type=int, default = 20, help='Min scWGS depth', required = False)
	parser.add_argument('--scDNACalls', type=str, default=1, help='scDNA calling file, "ground truth" (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--scDNAValidLong', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--scDNAValidSCom', type=str, default=1, help='scDNA support for SComatic calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--scDNA_supp_in_scRNA', type=str, default=1, help='SComatic support for scDNA calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--id', type=str, help='Sample ID', required = True)
	parser.add_argument('--outfile', type=str, help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	dp = args.dp
	scDNACalls = args.scDNACalls
	scDNAValid_L = args.scDNAValidLong
	scDNAValid_S = args.scDNAValidSCom
	RNASupp = args.scDNA_supp_in_scRNA
	sampleid = args.id
	out_prefix = args.outfile
	
	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Load df
	scDNACalls= pd.read_csv(scDNACalls, sep='\t', skiprows=29)
	scDNAValid_L = pd.read_csv(scDNAValid_L, sep='\t',  na_values=['.']).fillna(0)
	scDNAValid_S = pd.read_csv(scDNAValid_S, sep='\t',  na_values=['.']).fillna(0)
	RNASupp = pd.read_csv(RNASupp, sep='\t',  na_values=['.']).fillna(0)

	# 2. Create clinical annotation file
	Venn2(scDNAValid_L,scDNAValid_S,out_prefix)
	A,B,C,d = Venn3(dp, scDNAValid_L,scDNAValid_S,scDNACalls,out_prefix)
	scDNASupport(A,B,d,dp,scDNAValid_L,scDNAValid_S,sampleid,out_prefix)
	F1Score(A,B,C,d,RNASupp,sampleid,out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



