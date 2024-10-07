#!/usr/bin/python

import timeit
import argparse
import pandas as pd
import numpy as np
	
def variant_calling_step3(file,out_prefix,deltaVAF,deltaMCF,chrM_conta,min_ac_reads,min_ac_cells,clust_dist):

	#---
	# Command to focus only on high confidence cancer variants (HCCV)
	#---
	
	# Save comment lines
	f2 = open(out_prefix+ '.calling.step3.tsv','w')
	f3 = open(out_prefix+ '.calling.step3.unfiltered.tsv','w')
	with open(file, 'r') as f:
		for line in f:
			if line.startswith('#'):
				# Extract output file column names
				if '#CHROM' in line:
					line = line.rstrip('\n')
					output_column_names = line.split('\t')
				else:
					f2.write(line)
					f3.write(line)
			else:
				break
	f2.write('##INFO=FINAL_FILTER,Description=Final ilter status, including chrM contaminants and clustered sites\n')
	f3.write('##INFO=FINAL_FILTER,Description=Final ilter status, including chrM contaminants and clustered sites\n')	
	f2.close()
	f3.close()

	
	input_df = pd.read_csv(file, sep='\t',comment='#',names=output_column_names)
	input_df['INDEX'] = input_df['#CHROM'].astype(str) + ':' + input_df['Start'].astype(str) + ':' + input_df['ALT'].str.split(',', n=1, expand=True)[0]
	
	# Only keeping sites with alternative reads in cancer
	input_df = input_df[input_df['Cell_types']!='Non-Cancer']

	# Filtering multiallelic sites, keeping only one allele if it's 50x larger than the other
	# Otherwise deleting the line (it messes with filters later on)
	# This is important as SComatic will often detect low expr secondary alt allele in very high expr. (e.g. in chrM)
	input_df[['ALT', 'FILTER', 'Cell_types', 'Bc', 'Cc', 'VAF', 'MCF','MultiAllelic_filter']] =  input_df.apply(lambda x: MultiAllelic_filtering(x['REF'], x['ALT'], x['FILTER'], x['Cell_types'],x['Dp'], x['Nc'], x['Bc'], x['Cc'], x['VAF'], x['MCF'], x['Cancer'], x['Non-Cancer']), axis=1, result_type="expand") 
	input_df = input_df[input_df['MultiAllelic_filter']=='KEEP']
	input_df = input_df[output_column_names + ['INDEX']]

	#Special case for chrM due to contaminants:

	# Save chrM candidate SNVs to apply specific filters
	chrm_df = input_df[input_df['#CHROM']=='chrM'].copy()
	input_df = input_df[input_df['#CHROM']!='chrM']
	# Filtering homopolymeres, GnomAD/RNA editing/PON/clustered sites, 
	# multiallelic sites and sites with unsufficient # celltypes covered
	chrm_df = chrm_df[~chrm_df['FILTER'].str.contains('Min|LR|gnomAD|LC|RNA', regex=True)]
	if len(chrm_df)>0:
		# Applying deltaVAF and deltaMCF filters
		chrm_df['FILTER'] = chrm_df.apply(lambda x: 
			chrM_filtering(x['Cell_types'],x['Dp'],x['VAF'], x['MCF'],deltaVAF,deltaMCF), axis=1)


	# Filter 1: Non-cancer coverage filter
	input_df = input_df[~input_df['FILTER'].str.contains('Min_cell_types')]
	
	# Filter 2: Alt reads and cells in cancer filter
	input_df['Bc'] = [i.split(',')[1] if ',' in i else i for i in input_df['Bc']] #only keeping cancer info
	input_df['Cc'] = [i.split(',')[1] if ',' in i else i for i in input_df['Cc']] #only keeping cancer info
	input_df = input_df.astype({'Bc':'int','Cc':'int'})
	input_df = input_df[input_df['Bc']>=min_ac_reads]
	input_df = input_df[input_df['Cc']>=min_ac_cells]

	# Filter 3: Betabinomial filter in cancer cells
	input_df = input_df[~input_df['Cell_type_Filter'].str.contains(',Non-Significant|,Low-Significant', regex = True)]
	input_df = input_df[~input_df['Cell_type_Filter'].isin(['Non-Significant','Low-Significant'])]
	
	#Adding chrM SNVs detected above:
	unfiltered_df = pd.concat([input_df,chrm_df])
	unfiltered_df.to_csv(out_prefix+ '.calling.step3.unfiltered.tsv', sep='\t', index=False,  mode='a')

	# Filter 4: Noise filter:
	input_df = input_df[~input_df['FILTER'].str.contains('Noisy_site')] #multi allelic sites are filtered above
	
	# Filter 5: Homopolymer filter:
	input_df = input_df[~input_df['FILTER'].str.contains('LC_Upstream|LC_Downstream', regex = True)] 
	
	# Filter 6:RNA-Editing filter
	input_df = input_df[~input_df['FILTER'].str.contains('RNA_editing_db', regex = True)]

	# Filter 7: PoN filter
	input_df = input_df[~input_df['FILTER'].str.contains('PoN', regex = True)]
	
	# Filter 8: Betabinomial filter in non-cancer cells
	input_df = input_df[~input_df['Cell_type_Filter'].str.contains('Low-Significant,|PASS,', regex = True)]
	input_df = input_df[~input_df['FILTER'].str.contains('Cell_type_noise|Multiple_cell_types|Noisy_site', regex = True)]

	# Filter 9: gnomAD filter
	input_df = input_df[~input_df['FILTER'].str.contains('gnomAD', regex = True)]

	# Filter 10: Distance filter
	input_df['FILTER'] = tag_clustered_SNVs(input_df, clust_dist)
	input_df = input_df[~input_df['FILTER'].str.contains('dist', regex = True)]
	
	#Adding chrM SNVs detected above:
	input_df = pd.concat([input_df,chrm_df])

	#Filtering PASS SNVs (only left in chrM)
	filtered_df = input_df[input_df['FILTER'] =='PASS']

	# Write output
	filtered_df.to_csv(out_prefix+ '.calling.step3.tsv', sep='\t', index=False,  mode='a')
	
def chrM_filtering(CTYPES,DP,VAF,MCF,deltaVAFmin,deltaMCFmin):
	ctypes = CTYPES.split(',')
	if len(ctypes)>1:
		if ctypes[0]=='Cancer':
			i_Cancer = 0
			i_NonCancer = 1
		elif ctypes[1]=='Cancer':
			i_Cancer = 1
			i_NonCancer = 0
		DP1,DP2=DP.split(',')
		#Only look at higher depth loci
		if int(DP1)<100 or int(DP2)<100: 
			return 'LowDepth'
		else: 
			VAFs = VAF.split(',')
			MCFs = MCF.split(',')
			VAFCancer = float(VAFs[i_Cancer])
			VAFNonCancer = float(VAFs[i_NonCancer])
			deltaVAF = VAFCancer-VAFNonCancer
		
			MCFCancer = float(MCFs[i_Cancer])
			MCFNonCancer = float(MCFs[i_NonCancer])
			deltaMCF = MCFCancer-MCFNonCancer
			
			# mtSNVs are variants with high VAF in cancer and low VAF in non-cancer cells
			if deltaVAF < deltaVAFmin:
				return 'LowDeltaVAF'
		
			elif deltaMCF < deltaMCFmin:
				return 'LowDeltaMCF'
			else:
				return 'PASS'

	elif len(ctypes)==1:
		if ctypes[0]!="Cancer":
			return 'Non-Cancer'
		elif int(DP)<100:
			return 'LowDepth'
		elif float(VAF)<0.05:
			return 'LowVAF'
		elif float(MCF)<0.05:
			return 'LowMCF'
		else:
			return 'PASS'

def MultiAllelic_filtering(REF, ALT, FILTER, CTYPES, DP, NC, BC, CC, VAF, MCF, CancerInfo, NonCancerInfo):
	ref_dict = {'A':0, 'C':1, 'T':2, 'G':3}
	i_ref = ref_dict[REF]
	if 'Multi-allelic' in FILTER or '|' in ALT:
		ctypes = CTYPES.split(',')
		if len(ctypes)>1:
			if ctypes[0]=='Cancer':
				i_Cancer = 0
				i_NonCancer = 1
			elif ctypes[1]=='Cancer':
				i_Cancer = 1
				i_NonCancer = 0

			BCS = CancerInfo.split('|')[3].split(':')[:4]
			BCS = [int(i) for i in BCS]
			BCS[i_ref] = 0 # setting REF to 0 as we only consider ALT
			MAX = max(BCS)
			index = np.argmax(BCS)
			BCS[index] = 0 # removing max to select next "max"
			MAX2 = max(BCS)
			if not(MAX2/MAX<0.05): # one alt 20x larger than the other)
				return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'DELETE'

			ALT_Cancer = 'ACTG'[index]
			BC_Cancer = int(CancerInfo.split('|')[3].split(':')[index])
			CC_Cancer = int(CancerInfo.split('|')[2].split(':')[index])
			VAF_Cancer = round(BC_Cancer/int(DP.split(',')[i_Cancer]),4)
			MCF_Cancer = round(CC_Cancer/int(NC.split(',')[i_Cancer]),4)

			ALT_NonCancer = 'ACTG'[index]
			BC_NonCancer = int(NonCancerInfo.split('|')[3].split(':')[index])
			CC_NonCancer = int(NonCancerInfo.split('|')[2].split(':')[index])
			VAF_NonCancer = round(BC_NonCancer/int(DP.split(',')[i_NonCancer]),4)
			MCF_NonCancer = round(CC_NonCancer/int(NC.split(',')[i_NonCancer]),4)

			ALT = ','.join([ALT_NonCancer,ALT_Cancer])
			BC = ','.join([str(BC_NonCancer),str(BC_Cancer)])
			CC = ','.join([str(CC_NonCancer),str(CC_Cancer)])
			VAF = ','.join([str(VAF_NonCancer),str(VAF_Cancer)])
			MCF = ','.join([str(MCF_NonCancer),str(MCF_Cancer)])

			FILTER = FILTER.replace('Multi-allelic,','')
			FILTER = FILTER.replace(',Multi-allelic','')
			FILTER = FILTER.replace('Multi-allelic','')

			return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'KEEP'
		
		elif len(ctypes)==1:
			if ctypes[0]=='Cancer':
				BCS = CancerInfo.split('|')[3].split(':')[:4]
				BCS = [int(i) for i in BCS]
				BCS[i_ref] = 0 # setting REF to 0 as we only consider ALT
				MAX = max(BCS)
				index = np.argmax(BCS)
				BCS[index] = 0 # removing max to select next "max"
				MAX2 = max(BCS)
				if not(MAX2/MAX<0.05): # one alt 20x larger than the other)
					return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'DELETE'

				ALT = 'ACTG'[index]
				BC = int(CancerInfo.split('|')[3].split(':')[index])
				CC = int(CancerInfo.split('|')[2].split(':')[index])
				VAF = round(BC/int(DP),4)
				MCF= round(CC/int(NC),4)

				FILTER = FILTER.replace('Multi-allelic,','')
				FILTER = FILTER.replace(',Multi-allelic','')
				FILTER = FILTER.replace('Multi-allelic','')

				return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'KEEP'
			else:
				return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'DELETE'
	else:	
		return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'KEEP'
	
def tag_clustered_SNVs(df, clust_dist):
	df2 = df[df['FILTER']=='PASS']
	idx = df2['INDEX']
	a=[]
	for i in idx:
		chr,pos,base=i.split(':')
		a.append((chr,pos,base))
	b = sorted(a, key=lambda x: (x[0],x[1]))

	trash=[]

	for (chr,pos,base), (chr2,pos2,base2) in zip(b, b[1:]):
		if chr==chr2:
			if chr == 'chrM':
				continue
			if abs(int(pos)-int(pos2))<clust_dist:
				trash.append(':'.join([chr,pos,base]))
				trash.append(':'.join([chr2,pos2,base2]))
	trash = set(trash)
	updated_filter= df.apply(lambda x : 
							   modify_filter(x['INDEX'],x['FILTER'], clust_dist, trash),
							   axis=1)
	print(updated_filter)
	return updated_filter

def modify_filter(INDEX, FILTER, clust_dist, trash):
	clustered = 'Clust_dist{}'.format(str(clust_dist))
	if INDEX in trash:
		if FILTER == 'PASS':
			return clustered
		else:
			return FILTER + ',' + clustered
	else:
		return FILTER

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to perform the scRNA somatic variant calling')
	parser.add_argument('--infile', type=str, help='Input file with all samples merged in a single tsv', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--deltaVAF', type=float, default = 0.3, help='Delta VAF between cancer and non-cancer cells', required = True)
	parser.add_argument('--deltaMCF', type=float, default = 0.3, help='Delta MCF (cancer cell fraction) between cancer and non-cancer cells', required = True)
	parser.add_argument('--chrM_contaminant', type=str, default = 'True', help='Use this option if chrM contaminants are observed in non-cancer cells', required = False)
	parser.add_argument('--min_ac_reads', type=int, default = 2, help='Minimum ALT reads', required = False)
	parser.add_argument('--min_ac_cells', type=int, default = 3, help='Minimum mutated cells', required = False)
	parser.add_argument('--clust_dist', type=int, default = 5, help='Minimum distance required between two consecutive SNVs', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	infile = args.infile
	out_prefix = args.outfile
	deltaVAF = args.deltaVAF
	deltaMCF = args.deltaMCF
	chrM_conta = args.chrM_contaminant
	min_ac_reads = args.min_ac_reads
	min_ac_cells = args.min_ac_cells
	clust_dist = args.clust_dist

	# 1.1: Step 3: Filter only PASS mutations, not within clust_dist of each other, and apply specific chrM filters
	print ('\n- Variant calling step 3\n')
	variant_calling_step3(infile,out_prefix,deltaVAF,deltaMCF,chrM_conta,min_ac_reads,min_ac_cells,clust_dist)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


