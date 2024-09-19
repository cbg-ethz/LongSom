import timeit
import argparse
import glob
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def PlotCellTypeReannotation(umap,reannotation,binary,id,out_prefix):
	patient = {'B486':'P1','B497':'P2','B500':'P3'}
	colors = {'cancer_cancer':'#2a9df4',
		   	  'noncancer_cancer':'#ef9a9a',
			  'cancer_noncancer':'#ef5350',
			  'noncancer_noncancer':'#6d8cd4'
			  }

	reannotation = pd.read_csv(reannotation, sep='\t')
	d = {}
	d['cancer_cancer'] = reannotation[(reannotation['Automated_Cell_type']=='HGSOC') & (reannotation['Reannotated_cell_type']=='Cancer')]['Index']
	d['cancer_noncancer'] = reannotation[(reannotation['Automated_Cell_type']=='HGSOC') & (reannotation['Reannotated_cell_type']=='NonCancer')]['Index']
	d['noncancer_cancer'] = reannotation[(reannotation['Automated_Cell_type']=='Non-HGSOC') & (reannotation['Reannotated_cell_type']=='Cancer')]['Index']
	d['noncancer_noncancer'] = reannotation[(reannotation['Automated_Cell_type']=='Non-HGSOC') & (reannotation['Reannotated_cell_type']=='NonCancer')]['Index']

	d_rev = {i:k for k,v in d.items() for i in v}
	
	umap = pd.read_csv(umap)
	umap = umap.astype({'UMAP_1':'float','UMAP_2':'float'})
	umap['Color'] = umap.apply(lambda x: bc_to_colors(x['Barcode'],colors,d_rev),  axis=1)
	umap['Color'] = pd.Categorical(umap['Color'], ['#ef5350', '#6d8cd4', '#2a9df4', '#ef9a9a'])
	umap = umap.sort_values('Color')
	plt.figure(figsize = (4,4), dpi=600)
	ax = sns.scatterplot(data=umap,
				x='UMAP_1',
				y='UMAP_2',
				hue = 'Color',
				palette = {'#2a9df4':'#2a9df4',
		   	  			   '#ef9a9a':'#ef9a9a',
						   '#ef5350':'#ef5350',
						   '#6d8cd4':'#6d8cd4'
						   },
				s=13,
				legend=None
				)
	plt.tight_layout()
	ax.figure.savefig(out_prefix + '.UMAP.png', dpi=600, bbox_inches='tight')
	plt.close()

	# Plotting the mean of covered  SNVs mutated per cell
	binary = pd.read_csv(binary, sep='\t',index_col=0, na_values=[3])
	boxplot = {'Patient':[],'Barcode':[],'Cell Categ':[],'Total SNVs Mutated':[],'Mean Mutated Fraction':[],'Colors':[],}
	for k in d:
		for bc in d[k]:
			boxplot['Patient'].append(patient[id])
			boxplot['Barcode'].append(bc)
			boxplot['Cell Categ'].append(k)
			boxplot['Total SNVs Mutated'].append(np.sum(binary[bc]))
			boxplot['Mean Mutated Fraction'].append(np.mean(binary[bc]))
			boxplot['Colors'].append(colors[k])

	boxplot = pd.DataFrame(boxplot)
	boxplot.to_csv(out_prefix + '.boxplot.tsv', sep='\t')

	umap['ColorBurden'] = umap['Barcode'].map(dict(zip(boxplot['Barcode'],boxplot['Total SNVs Mutated'])))
	umap = umap[umap['Color'].isin(['#2a9df4','#ef9a9a'])]

	plt.figure(figsize = (4,4), dpi=600)
	cmap= sns.color_palette("viridis", as_cmap=True)
	ax = sns.scatterplot(data=umap,
				x='UMAP_1',
				y='UMAP_2',
				hue = 'ColorBurden',
				style='Color',
    			markers=['o','s'],
				s=10,
				palette = cmap,
				legend=None
				)
	sm = plt.cm.ScalarMappable(cmap= cmap, norm=None)
	cbar = plt.colorbar(sm, ax = plt.gca())
	cbar.set_label('Mutated Fraction')
	plt.tight_layout()
	ax.figure.savefig(out_prefix + '.UMAPBurden.png', dpi=600, bbox_inches='tight')
	plt.close()

	TP = (len(d['cancer_cancer'])/(len(d['cancer_cancer'])+len(d['noncancer_cancer'])))*100
	FN = (len(d['noncancer_cancer'])/(len(d['cancer_cancer'])+len(d['noncancer_cancer'])))*100
	FP = (len(d['cancer_noncancer'])/(len(d['cancer_noncancer'])+len(d['noncancer_noncancer'])))*100
	TN = (len(d['noncancer_noncancer'])/(len(d['cancer_noncancer'])+len(d['noncancer_noncancer'])))*100
	array = [[TP,FP],[FN,TN]]

	sns.set(font_scale=2)
	df_cm = pd.DataFrame(array, index = ['Cancer','NonCancer'],
					  columns = ['Cancer','NonCancer'])
	print(df_cm)
	plt.figure(figsize = (6,6), dpi=600)
	ax = sns.heatmap(df_cm, annot=True, cmap='RdPu', vmax =100, annot_kws={"fontsize":25}, fmt='.1f',cbar=True)
	ax.set_title('Patient {}'.format(patient[id]), y=1.35)
	ax.set_ylabel('Marker-based\nAnnotation', size = 25)
	ax.set_xlabel('Reannotation', size = 25)
	ax.xaxis.set_label_coords(0.5, 1.3)
	ax.yaxis.set_label_coords(-0.2, 0.5)
	plt.tick_params(axis='both', which='major', labelbottom = False, bottom=False, top = False, labeltop=True)
	plt.tight_layout()
	ax.figure.savefig(out_prefix + '.Heatmap.png', dpi=600, bbox_inches='tight')
	plt.close()


def bc_to_colors(bc,colors,d_rev):
	try:
		return colors[d_rev[bc]]
	except KeyError:
		return None

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to plot somatic and germline support of PoN-flaged SNVs')
	parser.add_argument('--umap', type=str, default=1, help='Gene expression-derived UMAP coordinates of the cells in the sample', required = True)
	parser.add_argument('--reannotation', type=str, default=1, help='Reannotation file (obtained by CellTypeReannotation.py)', required = True)
	parser.add_argument('--binary', type=str, default=1, help='Binary cell-mutation matrix (obtained by SingleCellGenotype.py)', required = True)
	parser.add_argument('--id', type=str, default=1, help='Sample id', required = True)
	parser.add_argument('--outfile', type=str, default='PlotPoNComparison.png', help='Plot', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	umap = args.umap
	reannotation = args.reannotation
	binary = args.binary
	id = args.id
	out_prefix = args.outfile

	# 1. Create clinical annotation file
	PlotCellTypeReannotation(umap,reannotation,binary,id,out_prefix)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



