import timeit
import argparse
import pandas as pd

def write_PoNComparison(YesPoN_L,NoPoN_L ,YesPoN_S,NoPoN_S, id, outfile):
	frac_L = NoPoN_L/(YesPoN_L+NoPoN_L)
	frac_S = NoPoN_S/(YesPoN_S+NoPoN_S)
	sample_to_patient = {'B486': 'P1', 'B497': 'P2', 'B500':'P3'}
	with open(outfile, 'w') as f:
		f.write("SampleID\tMethod\tYesPoN\tNoPoN\tPercentage\n")
		f.write("{}\t{}\t{}\t{}\t{}\n".format(sample_to_patient[id],'LongSom',YesPoN_L,NoPoN_L,frac_L))
		f.write("{}\t{}\t{}\t{}\t{}\n".format(sample_to_patient[id],'SComatic',YesPoN_S,NoPoN_S,frac_S))


def PoNComparison(YesPoN,NoPoN):
	YesPoN = pd.read_csv(YesPoN, sep='\t', skiprows=29)
	NoPoN = pd.read_csv(NoPoN, sep='\t', skiprows=29)
	#geno = pd.read_csv(geno, sep='\t',  na_values=['.']).fillna(0)
	
	NoPoN_SNVs = [i for i in list(NoPoN['INDEX']) if i not in list(YesPoN['INDEX'])]
	YesPoN_SNVs = list(YesPoN['INDEX'])

	return len(YesPoN_SNVs),len(NoPoN_SNVs)
	

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to evalutate the incidence of the LR PoN')
	parser.add_argument('--YesPoN_L', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--NoPoN_L', type=str, default=1, help='SComatic base calling file, without LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--YesPoN_S', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--NoPoN_S', type=str, default=1, help='SComatic base calling file, without LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample ID', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	NoPoN_L = args.NoPoN_L
	YesPoN_L = args.YesPoN_L
	NoPoN_S = args.NoPoN_S
	YesPoN_S = args.YesPoN_S
	id = args.id
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Create clinical annotation file
	YesPoN_L,NoPoN_L = PoNComparison(YesPoN_L,NoPoN_L)

	YesPoN_S,NoPoN_S = PoNComparison(YesPoN_S,NoPoN_S)

	write_PoNComparison(YesPoN_L,NoPoN_L ,YesPoN_S,NoPoN_S, id, outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



