#!/usr/bin/python

import timeit
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def scRNAFilter(scDNA_file,scRNA_file):
	
	#---
	# Command to focus only on high confidence cancer variants (HCCV)
	#---
	
	# Save comment lines and extract index
	INDEX={}
	f2 = open(outfile,'w')
	with open(scRNA_file, 'r') as f:
		for line in f:
			if line.startswith('#'):
				# Extract output file column names
				if '#CHROM' in line:
					line = line.rstrip('\n')
					output_column_names = line.split('\t')
				else:
					f2.write(line)
			else:
				line = line.rstrip('\n')
				index = line.split('\t')[-1]
				INDEX.add(index)
	f2.close()

	scDNA = pd.read_csv(scDNA_file,sep='\t')
	scDNA['INDEX'] = scDNA['#CHROM'] + ':' + scDNA['Start'].astype(str) + ':' + scDNA['ALT'].str.split(',', n=1, expand=True)[0]
	scDNA = scDNA[scDNA['INDEX'].isin(INDEX)]
	return scDNA

def AggregateClones(scDNA,outfile):
	Tum_Clones = [i for i in scDNA.columns if 'NonTum' not in i].unique()
	NonTum_Clones = [i for i in scDNA.columns if 'NonTum' in i ].unique()
	# Compute VAF in Tumor Clones from "raw" data
	scDNA[['DP_Tum','VAF_Tum']] = scDNA.apply(
		lambda x: getVAFandDP(x['ALT'],*[x[clone] for clone in Tum_Clones])
		)
	# Compute VAF in Non-Tumor Clones from "raw" data
	scDNA[['DP_NoNTum','VAF_NonTum']] = scDNA.apply(
		lambda x: getVAFandDP(x['ALT'],*[x[clone] for clone in NonTum_Clones])
		)
	cmap = plt.get_cmap('Reds', 100)
	cmap.set_bad('grey')
	scDNA['TumorColor'] = scDNA['VAF_Tum'].apply(lambda x: matplotlib.colors.rgb2hex(cmap(x)))
	scDNA['NonTumorColor'] = scDNA['VAF_NonTum'].apply(lambda x: matplotlib.colors.rgb2hex(cmap(x)))
	scDNA.to_csv(outfile, sep = '\t', index = False)

def getVAFandDP(ALT_base,*args):
	bases = {'A':0,'C':0,'T':0,'G':0}
	for clone in args:
		A,C,T,G = clone.split('|')[3].split(':')[:4]
		bases['A']+=int(A)
		bases['T']+=int(T)
		bases['C']+=int(C)
		bases['G']+=int(G)

		ALT = bases[ALT_base]
		DP = sum(bases.values())
		if ALT!=0:
			VAF = ALT/DP
		else:
			VAF = '.'
		return DP,VAF

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to compute scDNA clones VAF in scRNA-detected SNV candidates')
	parser.add_argument('--sc', type=str, help='Input file with all samples merged in a single tsv', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--editing', type=str, help='RNA editing file to be used to remove RNA-diting sites', required = False)
	parser.add_argument('--pon_SR', type=str, help='Short-read (SR) Panel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts', required = False)	
	parser.add_argument('--pon_LR', type=str, help='Long-read (LR) Panel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts', required = False)	
	parser.add_argument('--min_distance', type=int, default = 5, help='Minimum distance allowed between potential somatic variants [Default: 5]', required = False)
	parser.add_argument('--gnomAD_db', type=str, help='gnomAD v4 database file', required = False)
	parser.add_argument('--gnomAD_max', type=float, default = 0.01, help='Maximum gnomAD population VAF [default 0.01]', required = False)

	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	scDNA = args.scDNA
	scRNA = args.scRNA
	out_prefix = args.outfile

	# 1: Filter positions found in scRNA
	print ('\n- scDNA Clones Genotyping\n')
	
	scDNA = scRNAFilter(scDNA,scRNA)

	# 2: Aggregate clones
	print ('\n- scDNA Clones Aggregation\n')
	outfile = out_prefix + '.scDNAClonesGenotyping.tsv'
	AggregateClones(scDNA,outfile)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


