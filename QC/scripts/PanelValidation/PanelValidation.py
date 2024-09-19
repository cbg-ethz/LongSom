import timeit
import argparse
import pandas as pd
import glob

def PlotPanelSupport(panel_dir,longsom_dir,scomatic_dir,outfile):
	dfs = []
	for tsv in glob.glob(panel_dir+'/*CloneGenotype.tsv'):
		df = pd.read_csv(tsv, sep='\t', na_values=['.']).fillna(0)
		dfs.append(df)
	df_panel = pd.concat(dfs)

	dfs = []
	for tsv in glob.glob(longsom_dir+'/*.calling.step3.tsv'):
		df = pd.read_csv(tsv, sep='\t', skiprows=29)
		dfs.append(df)
	df_longsom = pd.concat(dfs)

	dfs = []
	for tsv in glob.glob(scomatic_dir+'/*.calling.step3.tsv'):
		df = pd.read_csv(tsv, sep='\t', skiprows=29)
		dfs.append(df)
	df_scomatic = pd.concat(dfs)

	df_panel['INDEX'] = df_panel['INDEX'].str.slice(0,-2)
	df_longsom['INDEX'] = df_longsom['INDEX'].str.slice(0,-2)
	df_scomatic['INDEX'] = df_scomatic['INDEX'].str.slice(0,-2)

	Panel_detected_longsom = [i for i in list(df_panel['INDEX']) if i in list(df_longsom['INDEX'])]
	Panel_detected_scomatic = [i for i in list(df_panel['INDEX']) if i in list(df_scomatic['INDEX'])]
	
	Panel_Germline_Support = list(df_panel[df_panel['NonCancer_VAF']>0]['INDEX'])
	Panel_Somatic_Support = [i for i in list(df_panel[df_panel['Cancer_VAF']>0]['INDEX']) if i not in Panel_Germline_Support]
	Panel_NoCov = list(df_panel[df_panel['Cancer_DP']==0]['INDEX'])
	Panel_LowCov = list(df_panel[df_panel['Cancer_DP'].between(1,5)]['INDEX'])

	with open (outfile,'w') as f:
		f.write('Panel_detected_longsom : {}\n'.format(len(Panel_detected_longsom)))
		for i in Panel_detected_longsom:
			f.write(f"{i}\n")

		f.write('Panel_detected_scomatic : {}\n'.format(len(Panel_detected_scomatic)))
		for i in Panel_detected_scomatic:
			f.write(f"{i}\n")

		f.write('Panel_Germline_Support : {}\n'.format(len(Panel_Germline_Support)))
		f.write('Panel_Somatic_Support : {}\n'.format(len(Panel_Somatic_Support)))
		f.write('Panel_LowCov : {}\n'.format(len(Panel_LowCov)))
		f.write('Panel_NoCov : {}\n'.format(len(Panel_NoCov)))
						  
def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get clinical and gene-level annotations')
	parser.add_argument('--panel_dir', type=str, default = 20, help='', required = False)
	parser.add_argument('--longsom_dir', type=str, default=1, help='', required = True)
	parser.add_argument('--scomatic_dir', type=str, default=1, help='', required = True)
	parser.add_argument('--outfile', type=str, default=1, help='', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	panel_dir = args.panel_dir
	longsom_dir = args.longsom_dir
	scomatic_dir = args.scomatic_dir
	outfile = args.outfile

	# 1.
	PlotPanelSupport(panel_dir,longsom_dir,scomatic_dir,outfile)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



