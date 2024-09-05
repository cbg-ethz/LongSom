
import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def PlotDistComparison(indir,outfile):
	dfs = []
	for tsv in glob.glob(indir+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t')
		dfs.append(df)

	df_plot = pd.concat(dfs)
	df_plot = pd.melt(df_plot, id_vars=['SampleID','DistFilter'], 
				   				value_vars=['FRAC_SUP','FRAC_GERM'], 
								var_name='scDNA Support',
								value_name='Fraction')
	df_plot['scDNA Support'] = df_plot['scDNA Support'].map({'FRAC_SUP':'Fraction of SNVs\nSupported',
														  	 'FRAC_GERM':'Fraction of SNV\nSupported as Germline',
														  })


	ax = sns.boxplot(
		x='scDNA Support',
		y='Fraction', 
		hue='DistFilter', 
		data=df_plot, 
		boxprops=dict(alpha=.5),
		showfliers=False,
	)

	sns.stripplot(
		x='scDNA Support',
		y='Fraction', 
		hue='DistFilter', 
		edgecolor='white',
		legend=None, 
		data=df_plot, 
		dodge=True, 
		alpha=0.8,
		s=10, 
		ax=ax
	)

	ax.set_ylabel("Fraction of SNV loci\nmutated",fontsize=18)
	ax.set_xlabel("scDNA Support",fontsize=18)
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
	PlotDistComparison(indir,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



