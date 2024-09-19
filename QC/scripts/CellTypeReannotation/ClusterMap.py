import argparse
import timeit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns



def CallClusterMap(bin, ctypes, height, width, out_prefix):
	
	ctypes = pd.read_csv(ctypes, sep = '\t')
	ctypes['CancerColorBeforeReannotation'] = ctypes['Automated_Cell_type'].apply(lambda x: '#ADBEE6' if x=='Non-HGSOC' else '#FFAEA5')
	data = pd.read_csv(bin, sep = '\t', index_col=0, na_values=[3]).fillna(0)
	colors_row = ['#067FD0' if '--' in i else '#E63B60' for i in data.index]
	
	ctypes = ctypes.sort_values(by=['Automated_Cell_type'], ascending = False)
	data = data[[i for i in ctypes['Index'] if i in data.columns]]
	out_file = out_prefix + '.ClusterMap.NoReannotation.pdf'
	colors_col = list(ctypes['CancerColorBeforeReannotation'])
	ClusterMap(data, height, width, colors_col, colors_row, out_file)
	
	ctypes = ctypes.sort_values(by=['Reannotated_cell_type'], ascending = False)
	data = data[[i for i in ctypes['Index'] if i in data.columns]]
	out_file = out_prefix + '.ClusterMap.Reannotation.pdf'
	colors_col = list(ctypes['Cancer_Color'])
	ClusterMap(data, height, width, colors_col, colors_row, out_file)

def ClusterMap(data, height, width, colors_col, colors_row, out_file):

	cmap = plt.get_cmap('Reds', 2)
	cmap.set_over('green')
	cmap.set_bad('white')
	cm = sns.clustermap(
		data,
		square=False,
		vmin=0, vmax=1,
		cmap=cmap, fmt='',
		linewidths=0, linecolor='lightgray',
		col_colors=colors_col,
		row_colors=colors_row,
		colors_ratio=[0.02,0.03],
		col_cluster=False,
		row_cluster=False,
		figsize=(width, height),
		xticklabels=False,
		yticklabels=False
	)

	cm.cax.set_visible(False)
	cm.ax_row_dendrogram.set_visible(False)
	cm.ax_heatmap.tick_params(right=False, bottom=False)
	cm.ax_heatmap.spines['top'].set_visible(True)
	cm.ax_heatmap.spines['right'].set_visible(True)
	cm.ax_heatmap.spines['bottom'].set_visible(True)
	cm.ax_heatmap.spines['left'].set_visible(True)
	# cm.gs.set_height_ratios([0.05, 0.05, 0.95])

	cm.savefig(out_file, dpi=100)
	plt.close()

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to generate clustermaps of cells x mutations')
	parser.add_argument('--bin', type=str, default=1, help='SComatic binary matrix (obtained by SingleCellGenotype.py)', required = True)
	parser.add_argument('--ctypes', type=str, help='Files linking cell-types to barcodes', required = True)
	parser.add_argument('--height', type=int, default=10, help='Out directory prefix', required = False)
	parser.add_argument('--width', type=int, default=15, help='Out directory prefix', required = False)
	parser.add_argument('--out_dir', type=str, help='Out directory prefix', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	bin = args.bin
	ctypes = args.ctypes
	height = args.height
	width = args.width
	out_prefix = args.out_dir

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Create clinical annotation file
	CallClusterMap(bin, ctypes, height, width, out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')
	