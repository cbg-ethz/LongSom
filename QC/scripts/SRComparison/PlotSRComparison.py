
import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def PlotSRComparison(indir,outfile):
	dfs = []
	for tsv in glob.glob(indir+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t')
		dfs.append(df)

	df_plot = pd.concat(dfs)
	d = {'LR':{},'SR': {}, 'Common':{}}
	lists = []
	for seq in ['LR','SR','Common']:
		for categ in ['PASS','GERM','SUPP','NOTUM_HIGHCOV']:
			d[seq][categ] = df_plot[df_plot['Sequencing Technology']==seq][categ].sum()
		TOT = sum(d[seq].values())
		FRAC_SUP = (d[seq]['PASS'])/TOT
		FRAC_GERM= (d[seq]['GERM'])/TOT

		lists.append([seq,'Fraction of SNVs\nSupported as Somatic',FRAC_SUP])
		lists.append([seq,'Fraction of SNVs\nSupported as Germline',FRAC_GERM])


	headers =['Sequencing Technology','scDNA Support','Fraction']
	
	df_plot = pd.DataFrame(lists, columns=headers)


	ax = sns.barplot(
		x='scDNA Support',
		y='Fraction', 
		hue='Sequencing Technology', 
		data=df_plot, 
		ci="sd", 
		edgecolor="black",
		errcolor="black",
		errwidth=1.5,
		capsize = 0.1,
		alpha=0.5
	)

	sns.stripplot(
		x='scDNA Support',
		y='Fraction', 
		hue='Sequencing Technology', 
		edgecolor='white',
		legend=None, 
		data=df_plot, 
		dodge=True, 
		alpha=0.8,
		s=10, 
		ax=ax
	)

	ax.set_ylabel("Fraction of scRNA SNVs\ndetected",fontsize=18)
	ax.set_xlabel("scWGS Support",fontsize=18)
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
	PlotSRComparison(indir,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



