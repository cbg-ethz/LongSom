#!/usr/bin/python

import timeit
import argparse
import pandas as pd
import numpy as np

def HCCV_SNV(SNVs,outfile,min_dp,deltaVAF,deltaMCF,clust_dist):

	#---
	# Command to focus only on high confidence cancer variants (HCCV)
	#---
	
	# Save comment lines
	f2 = open(outfile,'w')
	with open(SNVs, 'r') as f:
		for line in f:
			if line.startswith('#'):
				# Extract output file column names
				if '#CHROM' in line:
					line = line.rstrip('\n')
					output_column_names = line.split('\t')
				else:
					f2.write(line)
			else:
				break
	f2.write('##INFO=HCCV_FILTER,Description=Filter status of the variant site for cell reannotation (high-confidence cancer variants)\n')		
	f2.close()
	
	input_df = pd.read_csv(SNVs, sep='\t',comment='#',names=output_column_names)
	
	# INDEX CHR:POS:ALTBASE
	input_df['INDEX'] = input_df['#CHROM'].astype(str) + ':' + input_df['Start'].astype(str) + ':' + input_df['ALT'].str.split(',', n=1, expand=True)[0]
	
	# Only keeping sites with alternative reads in cancer
	input_df = input_df[input_df['Cell_types']!='Non-Cancer']

	# Filtering multiallelic sites, keeping only one allele if it's 50x larger than the other
	# Otherwise deleting the line (it messes with filters later on)
	# This is important as SComatic will often detect low expr secondary alt allele in very high expr. (e.g. in chrM) DP, NC, ALT, FILTER, CTYPES, BC, CC, VAF, MCF, CancerInfo, NonCancerInfo
	input_df[['ALT', 'FILTER', 'Cell_types', 'Bc', 'Cc', 'VAF', 'MCF','MultiAllelic_filter']] =  input_df.apply(lambda x: MultiAllelic_filtering(x['REF'], x['ALT'], x['FILTER'], x['Cell_types'],x['Dp'], x['Nc'], x['Bc'], x['Cc'], x['VAF'], x['MCF'], x['Cancer'], x['Non-Cancer']), axis=1, result_type="expand") 
	input_df = input_df[input_df['MultiAllelic_filter']=='KEEP']
	input_df = input_df[output_column_names + ['INDEX']]
	
	#Filter 1: DP filtering
	input_df['DP_FILTER'] = input_df.apply(lambda x: 
		DP_filtering(x['Cancer'],x['Non-Cancer'],min_dp), axis=1)
	input_df = input_df[input_df['DP_FILTER']=='PASS']
	input_df.to_csv(outfile+'2', sep='\t', index=False,  mode='a')

	#Special case for chrM due to contaminants:
	# Save chrM candidate SNVs to apply specific filters
	chrm_df = input_df[input_df['#CHROM']=='chrM'].copy()
	input_df = input_df[input_df['#CHROM']!='chrM']
	# Filtering homopolymeres, GnomAD/RNA editing/PON/clustered sites, 
	# multiallelic sites and sites with unsufficient # celltypes covered
	chrm_df = chrm_df[~chrm_df['FILTER'].str.contains('Min|LR|gnomAD|LC|RNA', regex=True)]

	# Filter 2: Noise filter:
	input_df = input_df[~input_df['FILTER'].str.contains('Noisy_site')] #multi allelic sites are filtered above
	
	# Filter 3: Homopolymer filter:
	input_df = input_df[~input_df['FILTER'].str.contains('LC_Upstream|LC_Downstream', regex = True)] 

	# Filter 4: gnomAD filter
	input_df = input_df[~input_df['FILTER'].str.contains('gnomAD', regex = True)]
	
	# Filter 5: RNA-Editing filter
	input_df = input_df[~input_df['FILTER'].str.contains('RNA_editing_db', regex = True)]

	# Filter 6: PoN filter
	input_df = input_df[~input_df['FILTER'].str.contains('PoN', regex = True)]

	#Adding chrM SNVs detected above:
	input_df = pd.concat([input_df,chrm_df])

	# Filter 7: PoN filterDelta VAF and MCF filtering
	input_df['HCCV_FILTER'] = input_df.apply(lambda x: 
		MCF_filtering(x['Cell_types'],x['VAF'], x['MCF'],deltaVAF,deltaMCF), axis=1)
	input_df.to_csv(outfile+'3', sep='\t', index=False,  mode='a')
	input_df = input_df[input_df['HCCV_FILTER']=='PASS']

	# Filter 8: Distance filter
	input_df['FILTER'] = tag_clustered_SNVs(input_df, clust_dist)
	input_df = input_df[~input_df['FILTER'].str.contains('dist', regex = True)]

	# Write output
	input_df.to_csv(outfile, sep='\t', index=False,  mode='a')

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
	idx = df['INDEX']
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


def DP_filtering(NonCancerInfo,CancerInfo,min_dp):
	try:
		DP1 = CancerInfo.split('|')[0]
		DP2 = NonCancerInfo.split('|')[0]
		if int(DP1)<min_dp or int(DP2)<min_dp:
			return 'LowDepth'
		else:
			return 'PASS'
	except AttributeError:
		return 'NoCov'

	
def MCF_filtering(CTYPES,VAF,MCF,deltaVAFmin,deltaMCFmin):
	ctypes = CTYPES.split(',')
	#If only 1 cell type is called
	if len(ctypes)==1 and ctypes[0]=="Cancer":
		if float(VAF) >= deltaVAFmin and float(MCF) >= deltaMCFmin:
			return 'PASS'
		else:
			return 'Low VAF/MCF'
	elif len(ctypes)>1:
		VAFs = VAF.split(',')
		MCFs = MCF.split(',')
		if ctypes[0]=='Cancer':
			VAFCancer = float(VAFs[0])
			VAFNonCancer = float(VAFs[1])
			MCFCancer = float(MCFs[0])
			MCFNonCancer = float(MCFs[1])
		elif ctypes[1]=='Cancer':
			VAFCancer = float(VAFs[1])
			VAFNonCancer = float(VAFs[0])
			MCFCancer = float(MCFs[1])
			MCFNonCancer = float(MCFs[0])
		
		if VAFCancer<0.05: 
			return 'NonSig'

		deltaVAF = VAFCancer-VAFNonCancer
		deltaMCF = MCFCancer-MCFNonCancer

		if VAFNonCancer>0.1 and deltaVAF<2*deltaVAFmin:
			return 'Heterozygous'
		
		if VAFNonCancer>0.2:
			return 'Heterozygous'

		# HCCV are variants with high VAF/MCF in cancer and low VAF/MCF in non-cancer cells
		# if deltaVAF < deltaVAFmin:
		#	return 'LowDeltaVAF'
		
		if deltaMCF < deltaMCFmin:
			return 'LowDeltaMCF'
		else:
			return 'PASS'
	else:
		return 'NonCancer'
	

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to perform High-Confidence Cancer Variants calling')
	parser.add_argument('--SNVs', type=str, help='', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--min_dp', type=float, default = 20, help='Minimum depth in both celltypes to call a HCCV', required = True)
	parser.add_argument('--deltaVAF', type=float, default = 0.1, help='Delta VAF between cancer and non-cancer cells', required = True)
	parser.add_argument('--deltaMCF', type=float, default = 0.4, help='Delta MCF (cancer cell fraction) between cancer and non-cancer cells', required = True)
	parser.add_argument('--clust_dist', type=int, default = 10000, help='Minimum distance required between two consecutive SNVs', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	SNVs = args.SNVs
	outfile = args.outfile
	min_dp = args.min_dp
	deltaVAF = args.deltaVAF
	deltaMCF = args.deltaMCF
	clust_dist = args.clust_dist

	# 1.2: Step 2: Add distance, editing and PoN filters
	print ('\n- High Confidence Cancer Variants calling\n')
	outfile = outfile + '.HCCV.tsv'
	HCCV_SNV(SNVs,outfile,min_dp,deltaVAF,deltaMCF,clust_dist)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


