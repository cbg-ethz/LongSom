import timeit
import argparse
import pandas as pd

def SRComparison(LR,SR,geno_LR,geno_SR,id,outfile):
	LR = list(pd.read_csv(LR, sep='\t', skiprows=29)['INDEX'])
	SR = list(pd.read_csv(SR, sep='\t', skiprows=29)['INDEX'])
	SR = ['chr'+str(i) if i!='MT' else 'chrM' for i in SR]
	L = [i for i in LR if i not in SR]
	S = [i for i in SR if i not in LR]
	C = [i for i in LR if i in SR]
	geno_LR = pd.read_csv(geno_LR, sep='\t',  na_values=['.']).fillna(0)
	geno_SR = pd.read_csv(geno_SR, sep='\t',  na_values=['.']).fillna(0)

	
	with open(outfile, 'w') as f:
		TITLE = ['SampleID','Sequencing Technology','#SNV','PASS','GERM','SUPP','NOTUM_HIGHCOV','FRAC_SUP','FRAC_GERM']
		f.write('\t'.join(TITLE)+'\n')
		write_data(f,id,geno_LR,L,'LR')
		write_data(f,id,geno_SR,S,'SR')
		write_data(f,id,geno_LR,C,'Common')



def genotyping(geno, SNVs):
	geno = geno[geno['INDEX'].isin(SNVs)]
	geno = geno[~geno['INDEX'].str.contains('chrM')]
	PASS = list(geno[(geno['Clone_Tum_MutatedStatus']=='PASS') & (geno['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX'])
	GERM = list(geno[geno['Clone_NonTum_MutatedStatus']=='PASS']['INDEX'])
	SUPP = [i for i in geno[(geno['Clone_Tum_VAF']>0) | geno['Clone_NonTum_VAF']>0]['INDEX'] if i not in PASS+GERM]
	NOTUM_HIGHCOV = [i for i in geno[geno['Clone_Tum_DP']>=17]['INDEX'] if i not in  PASS+GERM+SUPP]

	return len(PASS),len(GERM),len(SUPP),len(NOTUM_HIGHCOV)

def write_data(f,id,geno,SNVs,Condition):
	PASS,GERM,SUPP,NOTUM_HIGHCOV = genotyping(geno,SNVs)
	if sum([PASS,GERM,SUPP,NOTUM_HIGHCOV]) > 0:
		FRAC_SUP = (PASS+SUPP)/sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
		FRAC_GERM = (GERM)/sum([PASS,GERM,SUPP,NOTUM_HIGHCOV])
	else:
		FRAC_SUP = 0
		FRAC_GERM = 0
	LEN = len(SNVs)
	LIST = [id,Condition,LEN,PASS,GERM,SUPP,NOTUM_HIGHCOV,round(FRAC_SUP,2),round(FRAC_GERM,2)]
	f.write('\t'.join([str(i) for i in LIST])+'\n')

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--LR', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--SR', type=str, default=1, help='SComatic base calling file, without LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--geno_LR', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--geno_SR', type=str, default=1, help='scDNA support for LongSom calls (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample ID', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	SR = args.SR
	LR = args.LR
	geno_LR = args.geno_LR
	geno_SR = args.geno_SR
	id = args.id
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Compare
	SRComparison(LR,SR,geno_LR,geno_SR,id,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



