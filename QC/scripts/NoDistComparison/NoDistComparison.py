import timeit
import argparse
import pandas as pd

def DistComparison(YesDist,NoDist,genotype,id,outfile):
	YesDist = pd.read_csv(YesDist, sep='\t', skiprows=29)
	NoDist = pd.read_csv(NoDist, sep='\t', skiprows=29)
	genotype = pd.read_csv(genotype, sep='\t',  na_values=['.']).fillna(0)
	
	NoDist_SNVs = [i for i in NoDist['INDEX'] if i not in YesDist['INDEX']]
	YesDist_SNVs = list(YesDist['INDEX'])

	
	YesDist_germline_support_ratio = genotyping(genotype,YesDist_SNVs)
	NoDist_germline_support_ratio = genotyping(genotype,NoDist_SNVs)
	
	with open(outfile, 'w') as f:
		f.write("SampleID\tDist\tNo Dist\n")
		f.write("{}\t{}\t{}\n".format(id,YesDist_germline_support_ratio,NoDist_germline_support_ratio))
	
	
def genotyping(genotype, SNVs):
	genotype = genotype[genotype['INDEX'].isin(SNVs)]
	genotype_support = genotype[(genotype['Clone_Tum_MutatedStatus']=='PASS') & (genotype['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX']
	genotype_germline = genotype[genotype['Clone_NonTum_MutatedStatus']=='PASS']['INDEX']
	germline_support_ratio = len(genotype_germline)/len(genotype_support)
	return germline_support_ratio

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--YesDist', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--NoDist', type=str, default=1, help='SComatic base calling file, without LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--genotype', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample ID', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	NoDist = args.NoDist
	YesDist = args.YesDist
	genotype = args.genotype
	id = args.id
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Create clinical annotation file
	DistComparison(YesDist,NoDist,genotype,id,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



