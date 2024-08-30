import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def PlotMutationalBurden(indir,outfile):
	dfs = []
	for tsv in glob.glob(indir+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t', index_col=0)
		dfs.append(df)

	df_plot = pd.concat(dfs)
	df_plot.reset_index(drop=True, inplace=True)
	df_plot.sort_values('Patient', inplace = True)

	ax = sns.boxplot(
		x='Patient', 
		y='Mean Mutated Fraction', 
		hue='Colors', 
		data=df_plot, 
		boxprops=dict(alpha=.5),
		palette = {'#2a9df4':'#2a9df4',
					'#ef9a9a':'#ef9a9a',
					'#ef5350':'#ef5350',
					'#6d8cd4':'#6d8cd4'
					},
		showfliers=False,
		legend=None
	)

	sns.stripplot(
		x='Patient', 
		y='Mean Mutated Fraction', 
		hue='Colors',
		palette = {'#2a9df4':'#2a9df4',
					'#ef9a9a':'#ef9a9a',
					'#ef5350':'#ef5350',
					'#6d8cd4':'#6d8cd4'
					},
		legend=None, 
		data=df_plot, 
		dodge=True, 
		alpha=0.6, 
		ax=ax
	)
	ax.set_ylabel("Fraction of SNV loci\nmutated",fontsize=18)
	ax.set_xlabel("Patient",fontsize=18)
	plt.tight_layout()
	ax.figure.savefig(outfile, dpi=600, bbox_inches='tight',)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--indir', type=str, default=1, help='Directory where ', required = True)
	parser.add_argument('--outfile', type=str, default='PlotMutationalBurden.png', help='Plot', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	indir = args.indir
	outfile = args.outfile

	# 1. Create clinical annotation file
	PlotMutationalBurden(indir,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



