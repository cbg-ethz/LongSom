import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pywaffle import Waffle

def PlotscDNAValidation(indir,outfile):
	dfs = []
	for tsv in glob.glob(indir+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t', na_values=['.']).fillna(-1)
		dfs.append(df)

	df_plot = pd.concat(dfs)
	data = {}
	ALL = df_plot['INDEX']
	PASS = list(df_plot[(df_plot['Clone_Tum_MutatedStatus']=='PASS') & (df_plot['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX'])
	GERM = list(df_plot[df_plot['Clone_NonTum_MutatedStatus']=='PASS']['INDEX'])
	SUPP = [i for i in df_plot[(df_plot['Clone_Tum_VAF']>0) | df_plot['Clone_NonTum_VAF']>0]['INDEX'] if i not in PASS+GERM]
	NOTUM_LOWCOV = [i for i in df_plot[df_plot['Clone_Tum_DP'].between(1,9)]['INDEX'] if i not in  PASS+GERM+SUPP]
	NOTUM_HIGHCOV = [i for i in df_plot[df_plot['Clone_Tum_DP']>=10]['INDEX'] if i not in  PASS+GERM+SUPP]
	NOCOV = [i for i in df_plot[df_plot['Clone_Tum_DP']==0]['INDEX'] if i not in GERM]

	print(len(ALL))
	print(len(PASS))
	print(len(GERM))
	print(len(SUPP))
	print(len(NOTUM_LOWCOV ))
	print(len(NOTUM_HIGHCOV))
	print(len(NOCOV))
	print(len(ALL) - len(PASS+GERM+SUPP+NOTUM_LOWCOV+NOTUM_HIGHCOV+NOCOV))

	data['Somatic Call'] = len(PASS)
	data['Somatic Support'] = len(SUPP)
	data['Germline Call'] = len(GERM)
	data['No Support'] = len(NOTUM_HIGHCOV)
	data['No Support: Low Coverage'] = len(NOTUM_LOWCOV)
	data['No Coverage'] = len(NOCOV)

	data_div10 = {k:round(v/10) for k,v in data.items()}

	fig = plt.figure(
		FigureClass=Waffle,
		values=data_div10,
		colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', 'darkgrey', 'grey'],
		title={'label': 'scDNA Validation', 'loc': 'left'},
		labels=["{} ({:2.1%})".format(k,v/sum(data.values())) for k, v in data.items()],
		legend={'loc': 'lower left', 'bbox_to_anchor': (1, 0.5), 'ncol': 1, 'framealpha': 1},
		starting_location='NW',
		vertical=True,
		block_arranging_style='snake'
	)

	plt.tight_layout()
	plt.savefig(outfile,dpi=600)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--indir', type=str, default=1, help='Directory with PoN germline support tsv (obtained by PoNComparison.py)', required = True)
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
	PlotscDNAValidation(indir,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



