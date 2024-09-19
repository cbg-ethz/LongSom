import timeit
import argparse
import pandas as pd

def write_PoNComparison(PoN_SNVs_L,SNVs_L, PoN_SNVs_S, SNVs_S, id, outfile):
	Frac_L = PoN_SNVs_L/SNVs_L
	Frac_S = PoN_SNVs_S/SNVs_S
	sample_to_patient = {'B486': 'P1', 'B497': 'P2', 'B500':'P3'}
	with open(outfile, 'w') as f:
		f.write("SampleID\tMethod\tSNVs\tInPoN\tPercentage\n")
		f.write("{}\t{}\t{}\t{}\t{}\n".format(sample_to_patient[id],'LongSom',SNVs_L,PoN_SNVs_L,Frac_L))
		f.write("{}\t{}\t{}\t{}\t{}\n".format(sample_to_patient[id],'SComatic',SNVs_S,PoN_SNVs_S,Frac_S))


def PoNComparison(method,PoN):
	method = pd.read_csv(method, sep='\t', skiprows=29)
	PoN = pd.read_csv(PoN, sep='\t', skiprows=3)

	PoN['INDEX'] = PoN['#CHROM'].astype(str) + ':' + PoN['POS'].astype(str)
	method['INDEX'] = method['INDEX'].str.slice(0,-2)
	print(PoN['INDEX'][:5])
	print(method['INDEX'][:5])

	PoN_SNVs = [i for i in list(method['INDEX']) if i in list(PoN['INDEX'])]

	return len(PoN_SNVs),len(method)
	

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to evalutate the incidence of the LR PoN')
	parser.add_argument('--LongSom', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--SComatic', type=str, default=1, help='SComatic base calling file, with LR Panel of Normals (obtained by BaseCellCalling.step3.py)', required = True)
	parser.add_argument('--PoN', type=str, default=1, help='PoN file with all samples', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample ID', required = True)
	parser.add_argument('--outfile', help='Out prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	LongSom = args.LongSom
	SComatic = args.SComatic
	PoN = args.PoN
	id = args.id
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Create clinical annotation file
	PoN_SNVs_L,SNVs_L = PoNComparison(LongSom,PoN)

	PoN_SNVs_S,SNVs_S = PoNComparison(SComatic,PoN)

	write_PoNComparison(PoN_SNVs_L,SNVs_L, PoN_SNVs_S, SNVs_S, id, outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



