import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def PlotPoNComparison(indir):
	dfs = []
	for tsv in glob.glob(indir+'/*F1Scores.tsv'):
		df = pd.read_csv(tsv, sep='\t')
		dfs.append(df)

	df_plot = pd.concat(dfs)
	df_plot = pd.melt(df_plot, id_vars=['SampleID','Method'], 
				   				value_vars=['Precision','Sensitivity','F1'], 
								var_name='Metrics',
								value_name='Value')
	ax = sns.barplot(
		x='Metrics', 
		y='Value', 
		hue='Method', 
		data=df_plot, 
		ci="sd", 
		edgecolor="black",
		errcolor="black",
		errwidth=1.5,
		capsize = 0.1,
		alpha=0.5
	)
	sns.stripplot(
		x='Metrics', 
		y='Value', 
		hue='Method', 
		data=df_plot, dodge=True, alpha=0.6, ax=ax
	)
	# remove extra legend handles
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles[2:], labels[2:], title='Method', bbox_to_anchor=(1, 1.02), loc='upper left')
	plt.tight_layout()
	ax.figure.savefig(indir+'/F1_plot.png', dpi=600, bbox_inches='tight',)
	plt.close()

def PlotscDNASupport(indir):
	dfs = []
	for tsv in glob.glob(indir+'/*scDNASupport.tsv'):
		df = pd.read_csv(tsv, sep='\t')
		dfs.append(df)

	df_plot = pd.concat(dfs)
	df_plot = pd.melt(df_plot, id_vars=['SampleID','Method'], 
				   				value_vars=['FRAC_SUP','FRAC_GERM'], 
								var_name='scDNA Support',
								value_name='Fraction')
	df_plot['scDNA Support'] = df_plot['scDNA Support'].map({'FRAC_SUP':'Fraction of SNVs\nSupported',
														  	 'FRAC_GERM':'Fraction of SNV\nSupported as Germline',
														  })


	ax = sns.boxplot(
		x='scDNA Support',
		y='Fraction', 
		hue='Method', 
		data=df_plot, 
		boxprops=dict(alpha=.5),
		showfliers=False,
	)

	sns.stripplot(
		x='scDNA Support',
		y='Fraction', 
		hue='Method', 
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
	ax.figure.savefig(indir+'/scDNASupport.png', dpi=600, bbox_inches='tight',)
	plt.close()

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--indir', type=str, default=1, help='Directory with PoN grmline support tsv (obtained by PoNComparison.py)', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	indir = args.indir

	# Set outfile name
	print("Outfile: " , indir ,  "\n") 

	# 1. Create clinical annotation file
	PlotPoNComparison(indir)
	PlotscDNASupport(indir)
	
	

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



