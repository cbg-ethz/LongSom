import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from pywaffle import Waffle

def PlotWaffle(indir1,outfile1):
	dfs = []
	for tsv in glob.glob(indir1+'/*.tsv'):
		df = pd.read_csv(tsv, sep='\t', na_values=['.']).fillna(-1)
		dfs.append(df)

	df_plot = pd.concat(dfs)
	data = {}
	df_plot[df_plot['Clone_Tum_DP']>=10]
	ALL = df_plot['INDEX']
	PASS = list(df_plot[(df_plot['Clone_Tum_MutatedStatus']=='PASS') & (df_plot['Clone_NonTum_MutatedStatus']!='PASS')]['INDEX'])
	GERM = list(df_plot[df_plot['Clone_NonTum_MutatedStatus']=='PASS']['INDEX'])
	SOMSUPP = [i for i in df_plot[(df_plot['Clone_Tum_VAF']>0) & df_plot['Clone_NonTum_VAF']==0]['INDEX'] if i not in PASS+GERM]
	GERMSUPP = [i for i in df_plot[df_plot['Clone_NonTum_VAF']>0]['INDEX'] if i not in PASS+GERM]
	NOTUM_LOWCOV = [i for i in df_plot[df_plot['Clone_Tum_DP'].between(1,9)]['INDEX'] if i not in  PASS+GERM+SOMSUPP+GERMSUPP]
	NOTUM_HIGHCOV = [i for i in df_plot[df_plot['Clone_Tum_DP']>=10]['INDEX'] if i not in  PASS+GERM+SOMSUPP+GERMSUPP]
	NOCOV = [i for i in df_plot[df_plot['Clone_Tum_DP']==0]['INDEX'] if i not in GERM]

	print(len(ALL))
	print(len(PASS))
	print(len(GERM))
	print(len(SOMSUPP))
	print(len(GERMSUPP))
	print(len(NOTUM_LOWCOV ))
	print(len(NOTUM_HIGHCOV))
	print(len(NOCOV))
	print(len(ALL) - len(PASS+GERM+SOMSUPP+GERMSUPP+NOTUM_LOWCOV+NOTUM_HIGHCOV+NOCOV))

	data['Somatic Call'] = len(PASS)
	data['Somatic Support'] = len(SOMSUPP)
	data['Germline Call'] = len(GERM)
	data['Germline Support'] = len(GERMSUPP)
	data['No Support'] = len(NOTUM_HIGHCOV)

	#data_div10 = {k:round(v/10) for k,v in data.items()}

	fig = plt.figure(
		FigureClass=Waffle,
		values=data,
		columns=12,
		colors = ['#66c2a5', '#b5ead7', '#ff9aa2', '#ffb7b2', 'grey'],
		title={'label': 'scDNA Validation', 'loc': 'left'},
		labels=["{} ({:2.1%})".format(k,v/sum(data.values())) for k, v in data.items()],
		legend={'loc': 'lower left', 'bbox_to_anchor': (1, 0.4), 'ncol': 1, 'framealpha': 1},
		starting_location='NW',
		vertical=True,
		block_arranging_style='snake'
	)

	plt.tight_layout()
	plt.savefig(outfile1,dpi=600)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--indir1', type=str, default=1, help='Directory with PoN germline support tsv (obtained by PoNComparison.py)', required = True)
	parser.add_argument('--indir2', type=str, default=1, help='Directory with PoN germline support tsv (obtained by PoNComparison.py)', required = True)
	parser.add_argument('--outfile1', type=str, default='PlotPoNComparison.png', help='Plot', required = True)
	parser.add_argument('--outfile2', type=str, default='PlotPoNComparison.png', help='Plot', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	indir1 = args.indir1
	indir2 = args.indir2
	outfile1 = args.outfile1
	outfile2 = args.outfile2

	# Set outfile1 name
	print("Outfile: " , outfile1 ,  "\n") 

	# 1. Create clinical annotation file
	PlotWaffle(indir1,outfile1)
	PlotWaffle(indir2,outfile2)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



