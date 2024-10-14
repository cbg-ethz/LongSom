import timeit
import argparse
import pandas as pd
import numpy as np

def filter_input(bin,vaf,barcodes,min_cells_per_mut,min_pos_cov,out_prefix):
	bin = pd.read_csv(bin,sep='\t',index_col=0,na_values=[3,'.'])
	vaf = pd.read_csv(vaf,sep='\t',index_col=0,na_values=[3,'.'])
	barcodes = pd.read_csv(barcodes,sep='\t')
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
	
	# Filter all input files:
	vaf = vaf.loc[idx,cols]
	barcodes = barcodes[barcodes['Index'].isin(cols)]
	bin = pd.concat([bin,fusions_save[cols]])

	# Add reanno ctype colors
	barcodes['Cell_Reanno_Colors'] = barcodes['Reannotated_cell_type'].apply(lambda x: '#94C773' if x=='Non-Cancer' else '#8F79A1')
	
	# Write
	bin.to_csv(out_prefix + '.BinaryMatrix.tsv', sep='\t')
	vaf.to_csv(out_prefix + '.VAFMatrix.tsv', sep='\t')
	barcodes.to_csv(out_prefix + '.Barcodes.tsv', sep='\t', index = False)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to filter BnpC input matrix')
	parser.add_argument('--bin', type=str, default=1, help='SComatic binary matrix (obtained by SingleCellGenotype.py)', required = True)
	parser.add_argument('--vaf', type=str, default=1, help='SComatic VAF matrix (obtained by SingleCellGenotype.py)', required = True)
	parser.add_argument('--barcodes', type=str, default=1, help='Barcode to celltypes (obtained by CellTypeReannotation.py)', required = True)
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
	barcodes = args.barcodes
	min_cells_per_mut = args.min_cells_per_mut
	min_pos_cov = args.min_pos_cov
	out_prefix = args.outfile

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Create clinical annotation file
	filter_input(bin,vaf,barcodes,min_cells_per_mut,min_pos_cov,out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')