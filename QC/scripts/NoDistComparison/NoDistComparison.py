import timeit
import argparse
import pandas as pd

def DistComparison(YesDist,NoDist,geno,id,outfile):
	YesDist = pd.read_csv(YesDist, sep='\t', skiprows=29)
	NoDist = pd.read_csv(NoDist, sep='\t', skiprows=29)
	geno = pd.read_csv(geno, sep='\t',  na_values=['.']).fillna(0)
	
	if len(NoDist['INDEX']) - len(YesDist['INDEX']) <= 0:
		NoDist_SNVs = [i for i in list(NoDist['INDEX']) if i not in list(YesDist['INDEX'])]
		YesDist_SNVs = list(YesDist['INDEX'])
		for i in NoDist_SNVs:
			print(i)
	else:
		YesDist_SNVs = [i for i in list(YesDist['INDEX']) if i not in list(NoDist['INDEX'])]
		NoDist_SNVs = list(NoDist['INDEX'])
		for i in YesDist_SNVs:
			print(i)

	
	with open(outfile, 'w') as f:
		TITLE = ['SampleID','DistFilter','#DistFiltered','PASS','GERM','SUPP','NOTUM_HIGHCOV','FRAC_SUP','FRAC_GERM']
		f.write('\t'.join(TITLE)+'\n')
		write_data(f,id,geno,YesDist_SNVs,'Yes')
		write_data(f,id,geno,NoDist_SNVs,'No')


def genotyping(geno, SNVs):
	geno = geno[geno['INDEX'].isin(SNVs)]
	geno = geno[~geno['INDEX'].str.contains('chrM')]
	PASS = list(geno[(geno['Clone_Tum_MutatedStatus']=='PASS') & (geno['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX'])
	GERM = list(geno[geno['Clone_NonTum_MutatedStatus']=='PASS']['INDEX'])
	SUPP = [i for i in geno[(geno['Clone_Tum_VAF']>0) | geno['Clone_NonTum_VAF']>0]['INDEX'] if i not in PASS+GERM]
	NOTUM_HIGHCOV = [i for i in geno[geno['Clone_Tum_DP']>=10]['INDEX'] if i not in  PASS+GERM+SUPP]

	return len(PASS),len(GERM),len(SUPP),len(NOTUM_HIGHCOV)

def write_data(f,id,geno,Dist_SNVs,DistFilter):
	PASS,GERM,SUPP,NOTUM_HIGHCOV = genotyping(geno,Dist_SNVs)
	if sum([PASS,GERM,SUPP,NOTUM_HIGHCOV]) > 0:
		FRAC_SUP = (PASS+SUPP)/sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
		FRAC_GERM = (GERM)/sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
	else:
		FRAC_SUP = 0
		FRAC_GERM = 0
	LEN = len(Dist_SNVs)
	LIST = [id,DistFilter,LEN,PASS,GERM,SUPP,NOTUM_HIGHCOV,round(FRAC_SUP,2),round(FRAC_GERM,2)]
	f.write('\t'.join([str(i) for i in LIST])+'\n')

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--YesDist', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--NoDist', type=str, default=1, help='SComatic base calling file, without LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--geno', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample ID', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	NoDist = args.NoDist
	YesDist = args.YesDist
	geno = args.geno
	id = args.id
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Create clinical annotation file
	DistComparison(YesDist,NoDist,geno,id,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



