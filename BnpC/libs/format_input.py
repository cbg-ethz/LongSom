import timeit
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def filter_input(bin,vaf,scDNA,ctypes,min_cells_per_mut,min_pos_cov,out_prefix):
	bin = pd.read_csv(bin,sep='\t',index_col=0,na_values=[3,'.'])
	vaf = pd.read_csv(vaf,sep='\t',index_col=0,na_values=[3,'.'])
	scDNA = pd.read_csv(scDNA,sep='\t',na_values=['.'])
	ctypes = pd.read_csv(ctypes,sep='\t')
	#Save fusions:
	fusions = [i for i in bin.index if '--' in i]
	fusions_save = bin.loc[fusions,bin.columns]
	bin = bin.loc[[i for i in bin.index if i not in fusions]]
	
	#Count non-NaN non-zero values (mutated cells) per row
	bin = bin[bin.replace(0,np.nan).count(axis=1)>min_cells_per_mut]

	# Count cells with more than min_pos_cov positions covered
	bin = bin.loc[:,bin.count() > min_pos_cov]

	idx = list(bin.index) + fusions
	cols = bin.columns
	
	# Attributing a cmap color to scDNA VAF
	cmap = plt.get_cmap('Reds', 100)
	cmap.set_bad('grey')
	scDNA['TumorColor'] = scDNA['Clone_Tum_VAF'].apply(lambda x: matplotlib.colors.rgb2hex(cmap(x)))
	scDNA['NonTumorColor'] = scDNA['Clone_NonTum_VAF'].apply(lambda x: matplotlib.colors.rgb2hex(cmap(x)))
	for fusion in fusions:
		scDNA.loc[fusion]=pd.Series(dtype=float)
		scDNA.loc[fusion,'INDEX'] = fusion
		scDNA.loc[fusion,['TumorColor','NonTumorColor']] = 'blue'
	
	# Filter all input files:
	vaf = vaf.loc[idx,cols]
	if not 'chr' in idx[0]:
		scDNA['INDEX'] = [i.split('chr')[1] if '--' not in i else i for i in list(scDNA['INDEX'])]
	scDNA = scDNA[scDNA['INDEX'].isin(idx)]
	ctypes = ctypes[ctypes['Index'].isin(cols)]
	bin = pd.concat([bin,fusions_save[cols]])
	
	# Write
	bin.to_csv(out_prefix + '.BinaryMatrix.tsv', sep='\t')
	vaf.to_csv(out_prefix + '.VAFMatrix.tsv', sep='\t')
	scDNA.to_csv(out_prefix + '.scDNACloneGenotype.tsv', sep='\t', index = False)
	ctypes.to_csv(out_prefix + '.Barcodes.tsv', sep='\t', index = False)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to filter BnpC input matrix')
	parser.add_argument('--bin', type=str, default=1, help='SComatic binary matrix (obtained by SingleCellGenotype.py)', required = True)
	parser.add_argument('--vaf', type=str, default=1, help='SComatic VAF matrix (obtained by SingleCellGenotype.py)', required = True)
	parser.add_argument('--scDNA', type=str, default=1, help='SComatic scDNA clones VAF (obtained by scDNAClonesGenotyping.py)', required = True)
	parser.add_argument('--ctypes', type=str, default=1, help='Barcode to celltypes (obtained by CellTypeReannotation.py)', required = True)
	parser.add_argument('--min_cells_per_mut', type=int, default=5, help='SComatic+CellTypeReannotation base calling file (obtained by BaseCellCalling.step3.py)', required = False)
	parser.add_argument('--min_pos_cov', type=int, default=3, help='SComatic+CellTypeReannotation base calling file (obtained by BaseCellCalling.step3.py)', required = False)
	parser.add_argument('--outfile', default = 'Matrix.tsv', help='Out file', required = False)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	bin = args.bin
	vaf = args.vaf
	scDNA = args.scDNA
	ctypes = args.ctypes
	min_cells_per_mut = args.min_cells_per_mut
	min_pos_cov = args.min_pos_cov
	out_prefix = args.outfile

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Create clinical annotation file
	filter_input(bin,vaf,scDNA,ctypes,min_cells_per_mut,min_pos_cov,out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')