import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def PlotMutationalBurden(indir,outfile1,outfile2):
	dfs = []
	for tsv in glob.glob(indir+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t', index_col=0)
		dfs.append(df)

	df_plot = pd.concat(dfs)
	df_plot.reset_index(drop=True, inplace=True)
	df_plot['Cell Categ'] = pd.Categorical(df_plot['Cell Categ'], 
								['cancer_cancer', 'noncancer_cancer', 'noncancer_noncancer','cancer_noncancer'])
	df_plot.sort_values(['Patient','Cell Categ'], inplace = True)

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
	ax.figure.savefig(outfile1, dpi=600, bbox_inches='tight',)
	plt.close()

	ax = sns.boxplot(
		x='Patient', 
		y='Total SNVs Mutated', 
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
		y='Total SNVs Mutated', 
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
	ax.set_ylabel("SNV loci mutated",fontsize=18)
	ax.set_xlabel("Patient",fontsize=18)
	plt.tight_layout()
	ax.figure.savefig(outfile2, dpi=600, bbox_inches='tight',)
	plt.close()

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--indir', type=str, default=1, help='Directory where ', required = True)
	parser.add_argument('--outfile1', type=str, default='PlotMeanVAF.png', help='Plot', required = True)
	parser.add_argument('--outfile2', type=str, default='PlotMutationalBurden.png', help='Plot', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	indir = args.indir
	outfile1 = args.outfile1
	outfile2 = args.outfile2

	# 1. Create clinical annotation file
	PlotMutationalBurden(indir,outfile1,outfile2)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



