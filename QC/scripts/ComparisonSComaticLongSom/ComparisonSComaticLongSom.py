import timeit
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

def venn_diagram(SComatic,LongSom,scDNACalls,scDNAValidLong,scDNAValidSCom,scDNA_supp_in_scRNA,sampleid,out_prefix):
	SComatic = pd.read_csv(SComatic, sep='\t', skiprows=29)
	LongSom = pd.read_csv(LongSom, sep='\t', skiprows=29)

	# Defining mutation sets: without fusions, and only positions with suficient cov in RNA (scDNA cov is ensured by default)
	A = [i for i in SComatic['INDEX'] if '--' not in i]
	B = [i for i in LongSom['INDEX'] if '--' not in i]
	Ab = len([i for i in A if i not in B])
	aB = len([i for i in B if i not in A])
	AB = len([i for i in A if i in B])
	plt.figure(figsize=(4,4))
	v = venn3(subsets=(Ab,aB,AB),set_labels = ('SComatic','LongSom'))

	plt.savefig(out_prefix + '.Venn2.png', dpi=600)
	plt.close()


	scDNACalls= pd.read_csv(scDNACalls, sep='\t', skiprows=29)
	scDNAValidLong = pd.read_csv(scDNAValidLong, sep='\t',  na_values=['.']).fillna(0)
	scDNAValidSCom = pd.read_csv(scDNAValidSCom, sep='\t',  na_values=['.']).fillna(0)
	scDNA_supp_in_scRNA = pd.read_csv(scDNA_supp_in_scRNA, sep='\t',  na_values=['.']).fillna(0)

	scDNASupportedLongSom = scDNAValidLong[(scDNAValidLong['Clone_Tum_MutatedStatus']=='PASS') & (scDNAValidLong['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX']
	scDNASupportedGermlineLongSom = scDNAValidLong[scDNAValidLong['Clone_NonTum_MutatedStatus']=='PASS']['INDEX']
	LongSomCoveredInscDNA = scDNAValidLong[scDNAValidLong['Clone_Tum_DP']>5]['INDEX']

	scDNASupportedSComatic = scDNAValidSCom[(scDNAValidSCom['Clone_Tum_MutatedStatus']=='PASS') & (scDNAValidSCom['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX']
	scDNASupportedGermlineSComatic = scDNAValidSCom[scDNAValidSCom['Clone_NonTum_MutatedStatus']=='PASS']['INDEX']
	SComaticCoveredInscDNA = scDNAValidSCom[scDNAValidSCom['Clone_Tum_DP']>5]['INDEX']

	# Defining mutation sets: without fusions, and only positions with suficient cov in RNA (scDNA cov is ensured by default)
	A = [i for i in SComatic['INDEX'] if '--' not in i if i in list(SComaticCoveredInscDNA)]
	B = [i for i in LongSom['INDEX'] if '--' not in i if i in list(LongSomCoveredInscDNA)]
	C = list(scDNACalls['INDEX']) + list(scDNASupportedLongSom) + list(scDNASupportedSComatic)

	Abc = len([i for i in A if i not in set(B+C)])
	aBc = len([i for i in B if i not in set(A+C)])
	abC = len([i for i in C if i not in set(A+B)])
	ABc = len([i for i in A if i in set(B) and i not in set(C)])
	aBC = len([i for i in C if i in set(B) and i not in set(A)])
	AbC = len([i for i in A if i in set(C) and i not in set(B)])
	ABC = len([i for i in A if i in set(C) and i in set(B)])

	plt.figure(figsize=(4,4))
	v = venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),set_labels = ('SComatic','LongSom','scDNA'))

	plt.savefig(out_prefix + '.Venn3.png', dpi=600)

	# Defining Validation Set: SNVs detected in scDNA with at least 1 mutated read in scRNA 
	scDNAValidSet = list(scDNA_supp_in_scRNA[(scDNA_supp_in_scRNA['HGSOC_VAF']>0) | (scDNA_supp_in_scRNA['Non-HGSOC_VAF']>0)]['INDEX'])
	with open(out_prefix + '.ValidationSet.tsv','w') as f:
		f.write('ValidationSet\n')
		for i in scDNAValidSet:
			f.write(i + '\n')

	with open(out_prefix + '.F1Scores.tsv','w') as f:
		f.write('SampleID\tMethod\tPrecision\tSensitivity\tF1\n')

		### LongSom Precision:
		# True positives re. precision are SNVs detected by LongSom and detected in scDNA 
		# or with at least 1 supporting read in scDNA and no reads supporting any other alternative allele
		TP_prec = len([i for i in B if i in set(scDNAValidSet + list(scDNASupportedLongSom))])
		# False positives are SNVs detected in scDNA but not in LongSom, and SNVs detected in LongSom but germlines in scDNA
		FP = len( 
				set( [i for i in B if i not in set(scDNAValidSet)] + list(scDNASupportedGermlineLongSom) ) 
			)
		Precision = TP_prec / (TP_prec + FP)

		# LongSom Sensitivity:
		# True positives re. sensitivity are SNVs detected by LongSom 
		TP_sens = len([i for i in B if i in set(scDNAValidSet)])
		# False negatives are dtected in scDNA, but not by Longsom
		FN = len([i for i in set(scDNAValidSet) if i not in B])
		Sensitivity = TP_sens / (TP_sens + FN)

		# LongSom F1 score:
		F1 = (2*TP_prec) / ((2*TP_prec) + FP + FN) 
		
		f.write(f'{sampleid}\tLongSom\t{Precision}\t{Sensitivity}\t{F1}\n')

		### SComatic Precision:
		TP_prec = len([i for i in A if i in set(scDNAValidSet + list(scDNASupportedSComatic))])
		FP = len( 
				set( [i for i in A if i not in set(scDNAValidSet)] + list(scDNASupportedGermlineSComatic) ) 
			)
		Precision = TP_prec / (TP_prec + FP)

		# LongSom Sensitivity:
		TP_sens = len([i for i in A if i in set(scDNAValidSet)])
		FN = len([i for i in set(scDNAValidSet) if i not in A])
		Sensitivity = TP_sens / (TP_sens + FN)

		# SComatic F1 score:
		F1 = (2*TP_prec) / ((2*TP_prec) + FP + FN) 

		f.write(f'{sampleid}\tSComatic\t{Precision}\t{Sensitivity}\t{F1}\n')


def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--SComatic', type=str, default=1, help='SComatic base calling file (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--LongSom', type=str, default=1, help='SComatic+CellTypeReannotation base calling file (obtained by BaseCellCalling.step3.py)', required = True)
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

	SComatic = args.SComatic
	LongSom = args.LongSom
	scDNACalls = args.scDNACalls
	scDNAValidLong = args.scDNAValidLong
	scDNAValidSCom = args.scDNAValidSCom
	scDNA_supp_in_scRNA = args.scDNA_supp_in_scRNA
	sampleid = args.id
	out_prefix = args.outfile

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Create clinical annotation file
	venn_diagram(SComatic,LongSom,scDNACalls,scDNAValidLong,scDNAValidSCom,scDNA_supp_in_scRNA,sampleid,out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



