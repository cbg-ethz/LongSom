import timeit
import argparse
import pandas as pd

def PoNComparison(YesPoN,NoPoN,genotype,id,outfile):
	YesPoN = pd.read_csv(YesPoN, sep='\t', skiprows=29)
	NoPoN = pd.read_csv(NoPoN, sep='\t', skiprows=29)
	genotype = pd.read_csv(genotype, sep='\t',  na_values=['.']).fillna(0)
	
	NoPoN_SNVs = [i for i in NoPoN['INDEX'] if i not in YesPoN['INDEX']]
	YesPoN_SNVs = list(YesPoN['INDEX'])
	print('len(PoN_SNVs) :',len(YesPoN_SNVs))
	print('len(NoPoN_SNVs) :',len(NoPoN_SNVs))
	
	YesPoN_germline_tot_ratio,YesPoN_germline_support_ratio,supp,germ = genotyping(genotype,YesPoN_SNVs)
	print('YesPoN_germline_support_ratio : ',YesPoN_germline_support_ratio)
	print('PoN scDNA sup ', supp)
	print('PoN scDNA germ ', germ)
	NoPoN_germline_tot_ratio,NoPoN_germline_support_ratio,supp,germ = genotyping(genotype,NoPoN_SNVs)
	print('NoPoN_germline_support_ratio : ', NoPoN_germline_support_ratio)
	print('NoPoN scDNA sup ', supp)
	print('NoPoN scDNA germ ', germ)
	with open(outfile, 'w') as f:
		f.write("SampleID\tYes\tNo\t#PoNFiltered\n")
		f.write("{}\t{}\t{}\t{}\n".format(id,YesPoN_germline_tot_ratio,NoPoN_germline_tot_ratio,len(NoPoN_SNVs)))
	
def genotyping(genotype, SNVs):
	genotype = genotype[genotype['INDEX'].isin(SNVs)]
	genotype = genotype[(genotype['Clone_Tum_VAF']>0) | (genotype['Clone_NonTum_VAF']>0)]
	tot = len(genotype)
	genotype_support = genotype[(genotype['Clone_Tum_VAF']>0) & (genotype['Clone_NonTum_VAF']==0)]['INDEX']
	genotype_germline = genotype[genotype['Clone_NonTum_VAF']>0]['INDEX']
	if tot>0:
		germline_tot_ratio = len(genotype_germline)/tot
	else:
		germline_tot_ratio = 0
	if len(genotype_support):
		germline_support_ratio = len(genotype_germline)/len(genotype_support)
	else:
		germline_support_ratio = 1
	
	return germline_tot_ratio,germline_support_ratio,len(genotype_germline),len(genotype_support)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--YesPoN', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--NoPoN', type=str, default=1, help='SComatic base calling file, without LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--genotype', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample ID', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	NoPoN = args.NoPoN
	YesPoN = args.YesPoN
	genotype = args.genotype
	id = args.id
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Create clinical annotation file
	PoNComparison(YesPoN,NoPoN,genotype,id,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



