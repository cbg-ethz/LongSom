
import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def PlotPoNComparison(indir,outfile):
	dfs = []
	for tsv in glob.glob(indir+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t')
		dfs.append(df)

	df_plot = pd.concat(dfs)
	df_plot['Percentage'].round(2)
	df_plot = df_plot.sort_values(['SampleID'])
	df_plot = df_plot[df_plot['SampleID']!='P2']

	ax = sns.barplot(
		x='SampleID',
		y='Percentage', 
		hue='Method',
		data=df_plot, 
		ci="sd", 
		edgecolor="black",
		errcolor="black",
		errwidth=1.5,
		capsize = 0.1,
		alpha=0.5
	)

	ax.set_xlabel("Patient",fontsize=18)
	ax.set_ylabel("Fraction of SNV loci\nfiltered by matched PoN_LR",fontsize=18)

	# remove extra legend handles
	# handles, labels = ax.get_legend_handles_labels()
	# ax.legend(handles[2:], labels[2:], title='With PoN LR', bbox_to_anchor=(1, 1.02), loc='upper left')
	plt.tight_layout()
	ax.figure.savefig(outfile, dpi=600, bbox_inches='tight',)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--indir', type=str, default=1, help='Directory with PoN grmline support tsv (obtained by PoNComparison.py)', required = True)
	parser.add_argument('--outfile', type=str, default='PlotPoNComparison.png', help='Plot', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	indir = args.indir
	outfile = args.outfile

	# Set outfile name
	print("Outfile: " , outfile ,  "\n") 

	# 1. Create clinical annotation file
	PlotPoNComparison(indir,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



