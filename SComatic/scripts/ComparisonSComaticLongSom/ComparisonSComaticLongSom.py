import timeit
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

def venn_diagram(SComatic,LongSom,scDNACalls,scDNAValidLong,scDNAValidSCom,scDNA_supp_in_scRNA,out_prefix):
	SComatic = pd.read_csv(SComatic, sep='\t', skiprows=29)
	LongSom = pd.read_csv(LongSom, sep='\t', skiprows=29)
	scDNACalls= pd.read_csv(scDNACalls, sep='\t', skiprows=29)
	scDNAValidLong = pd.read_csv(scDNAValidLong, sep='\t',  na_values=['.']).fillna(-1)
	scDNAValidSCom = pd.read_csv(scDNAValidSCom, sep='\t',  na_values=['.']).fillna(-1)
	scDNA_supp_in_scRNA = pd.read_csv(scDNA_supp_in_scRNA, sep='\t',  na_values=['.']).fillna(-1)

	scDNASupportedLongSom = scDNAValidLong[(scDNAValidLong['Clone_Tum_VAF']>0) & (scDNAValidLong['Clone_NonTum_VAF']<=0)]['INDEX']
	scDNASupportedGermlineLongSom = scDNAValidLong[scDNAValidLong['Clone_NonTum_VAF']>0]['INDEX']

	scDNASupportedSComatic = scDNAValidSCom[(scDNAValidSCom['Clone_Tum_VAF']>0) & (scDNAValidSCom['Clone_NonTum_VAF']<=0)]['INDEX']
	scDNASupportedGermlineSComatic = scDNAValidSCom[scDNAValidSCom['Clone_NonTum_VAF']>0]['INDEX']

	A = list(SComatic['INDEX'])
	B = list(LongSom['INDEX'])
	C = list(scDNACalls['INDEX']) + list(scDNASupportedLongSom) + list(scDNASupportedSComatic)

	Abc = len([i for i in A if i not in set(B+C)])
	aBc =  len([i for i in B if i not in set(A+C)])
	abC =  len([i for i in C if i not in set(A+B)])
	ABc = len([i for i in A if i in set(B) and i not in set(C)])
	aBC = len([i for i in C if i in set(B) and i not in set(A)])
	AbC = len([i for i in A if i in set(C) and i not in set(B)])
	ABC = len([i for i in A if i in set(C) and i in set(B)])

	plt.figure(figsize=(4,4))
	v = venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),set_labels = ('SComatic','LongSom','scDNA'))

	plt.savefig(out_prefix + '.Venn3.png', dpi=600)

	# Defining Validation Set:
	scDNAValidSet = list(scDNA_supp_in_scRNA[(scDNA_supp_in_scRNA['HGSOC_VAF']>0) | (scDNA_supp_in_scRNA['Non-HGSOC_VAF']>0)]['INDEX'])
	print(len(scDNAValidSet))

	# LongSom Precision:
	TP_prec = len([i for i in B if i in set(scDNAValidSet + list(scDNASupportedLongSom))])
	FP = len( 
			set( [i for i in B if i not in set(scDNAValidSet)] + list(scDNASupportedGermlineLongSom) ) 
		  )
	Precision = TP_prec / (TP_prec + FP)

	# LongSom Sensitivity:
	TP_sens = len([i for i in B if i in set(scDNAValidSet)])
	FN = len([i for i in set(scDNAValidSet) if i not in B])
	Sensitivity = TP_sens / (TP_sens + FN)

	# LongSom F1 score:
	F1 = (2*TP_prec) / ((2*TP_prec) + FP + FN) 
	
	with open(out_prefix + '.LongSomScores.tsv','w') as f:
		f.write('Precision\tSensitivity\tF1\n')
		f.write(f'{Precision}\t{Sensitivity}\t{F1}\n')
	
	# SComatic Precision:
	TP_prec = len([i for i in A if i in set(scDNAValidSet + list(scDNASupportedSComatic))])
	FP = len( 
			set( [i for i in A if i not in set(scDNAValidSet)] + list(scDNASupportedGermlineSComatic) ) 
		  )
	Precision = TP_prec / (TP_prec + FP)

	# LongSom Sensitivity:
	TP_sens = len([i for i in A if i in set(scDNAValidSet)])
	FN = len([i for i in set(scDNAValidSet) if i not in A])
	Sensitivity = TP_sens / (TP_sens + FN)

	# SDComatic F1 score:
	F1 = (2*TP_prec) / ((2*TP_prec) + FP + FN) 
	
	with open(out_prefix + '.SComaticScores.tsv','w') as f:
		f.write('Precision\tSensitivity\tF1\n')
		f.write(f'{Precision}\t{Sensitivity}\t{F1}\n')


def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--SComatic', type=str, default=1, help='SComatic base calling file (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--LongSom', type=str, default=1, help='SComatic+CellTypeReannotation base calling file (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--scDNACalls', type=str, default=1, help='scDNA calling file, "ground truth" (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--scDNAValidLong', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--scDNAValidSCom', type=str, default=1, help='scDNA support for SComatic calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--scDNA_supp_in_scRNA', type=str, default=1, help='SComatic support for scDNA calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
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
	out_prefix = args.outfile

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Create clinical annotation file
	venn_diagram(SComatic,LongSom,scDNACalls,scDNAValidLong,scDNAValidSCom,scDNA_supp_in_scRNA,out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



