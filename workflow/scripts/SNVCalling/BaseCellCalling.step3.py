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

	# Only keeping sites called in cancer
	input_df = input_df[input_df['Cell_types']!='Non-Cancer']

	# Filtering multiallelic sites, keeping only one allele if it's 50x larger than the other
	# Otherwise deleting the line (it messes with filters later on)
	# This is important as SComatic will often detect low expr secondary alt allele in very high expr. (e.g. in chrM)
	input_df[['ALT', 'FILTER', 'Cell_types', 'Bc', 'Cc', 'VAF', 'MCF','STEP3FILTER']] =  input_df.apply(lambda x: MultiAllelic_filtering(x['REF'], x['ALT'], x['FILTER'], x['Cell_types'],x['Dp'], x['Nc'], x['Bc'], x['Cc'], x['VAF'], x['MCF'], x['Cancer'], x['Non-Cancer']), axis=1, result_type="expand") 
	#input_df = input_df[input_df['MultiAllelic_filter']=='KEEP']

	input_df['INDEX'] = input_df['#CHROM'].astype(str) + ':' + input_df['Start'].astype(str) + ':' + input_df['ALT'].str.split(',', n=1, expand=True)[0]

	#Special case for chrM due to contaminants:
	# Save chrM candidate SNVs to apply specific filters
	chrm_df = input_df[input_df['#CHROM']=='chrM'].copy()
	input_df = input_df[input_df['#CHROM']!='chrM']
	# Filtering homopolymeres, GnomAD/RNA editing/PON/clustered sites, 
	# multiallelic sites and sites with unsufficient # celltypes covered
	chrm_df = chrm_df[~chrm_df['FILTER'].str.contains('Min|LR|gnomAD|LC|RNA', regex=True)]
	if len(chrm_df)>0:
		# Applying deltaVAF and deltaMCF filters
		chrm_df['STEP3FILTER'] = chrm_df.apply(lambda x: 
			chrM_filtering(x['STEP3FILTER'],x['Cell_types'],x['Dp'],x['VAF'], x['MCF'],deltaVAF,deltaMCF), axis=1)

	# Filter 1: Non-cancer coverage filter
	input_df = input_df[~input_df['FILTER'].str.contains('Min_cell_types')]

	#Filter 2: Min. alt reads and cells in cancer
	input_df['STEP3FILTER'] = input_df.apply(lambda x: BC_CC_filtering(x['STEP3FILTER'],x['ALT'],x['Cancer'],min_ac_reads,min_ac_cells), axis=1)

	# Filter 3 and 8: Betabinomial filtering in cancer cells (Filter 3) and non-cancer cells (Filter 8)
	input_df['STEP3FILTER'] = input_df.apply(lambda x: BetaBino_filtering(x['STEP3FILTER'],x['Cell_types'],x['Cell_type_Filter'],x['Start']), axis=1)

	# Filter 4: Noise filter:
	input_df = input_df[~input_df['FILTER'].str.contains('Noisy_site')] #multi allelic sites are filtered above
	
	# Filter 5: Homopolymer filter:
	input_df = input_df[~input_df['FILTER'].str.contains('LC_Upstream|LC_Downstream', regex = True)] 
	
	# Filter 6:RNA-Editing filter
	input_df = input_df[~input_df['FILTER'].str.contains('RNA_editing_db', regex = True)]

	# Filter 7: PoN filter
	input_df = input_df[~input_df['FILTER'].str.contains('PoN', regex = True)]
	
	# Filter 8: Betabinomial filter in non-cancer cells
	input_df = input_df[~input_df['FILTER'].str.contains('Cell_type_noise', regex = True)]

	# Filter 9: gnomAD filter
	input_df = input_df[~input_df['FILTER'].str.contains('gnomAD', regex = True)]

	#Adding chrM SNVs detected above:
	input_df = pd.concat([input_df,chrm_df])

	# Filter 10: Distance filter
	input_df['STEP3FILTER'] = tag_clustered_SNVs(input_df, clust_dist)
	input_df.to_csv(out_prefix+ '.calling.step3.unfiltered.tsv', sep='\t', index=False,  mode='a')
	input_df = input_df[~input_df['STEP3FILTER'].str.contains('dist', regex = True)]

	#Filtering PASS SNVs (only left in chrM)
	filtered_df = input_df[input_df['STEP3FILTER'] =='PASS']

	# Write output
	filtered_df.to_csv(out_prefix+ '.calling.step3.tsv', sep='\t', index=False,  mode='a')

	
def chrM_filtering(STEP3FILTER,CTYPES,DP,VAF,MCF,deltaVAFmin,deltaMCFmin):
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
			if STEP3FILTER == 'PASS':
				return 'LowDepth'
			else:
				return STEP3FILTER+',LowDepth'
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
				if STEP3FILTER == 'PASS':
					return 'LowDeltaVAF'
				else:
					return STEP3FILTER+',LowDeltaVAF'

		
			elif deltaMCF < deltaMCFmin:
				if STEP3FILTER == 'PASS':
					return 'LowDeltaMCF'
				else:
					return STEP3FILTER+',LowDeltaMCF'
			else:
				return STEP3FILTER

	elif len(ctypes)==1:
		if int(DP)<100:
			if STEP3FILTER == 'PASS':
				return 'LowDepth'
			else:
				return STEP3FILTER+',LowDepth'
		elif float(VAF)<0.05:
			if STEP3FILTER == 'PASS':
				return 'LowVAF'
			else:
				return STEP3FILTER+',LowVAF'
		elif float(MCF)<0.05:
			if STEP3FILTER == 'PASS':
				return 'LowMCF'
			else:
				return STEP3FILTER+',LowMCF'
		else:
			return STEP3FILTER

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

			if not(MAX2/MAX<0.05): # one alt 20x larger than the other)
				return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'Multi-Allelic'

			return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'PASS'
		
		elif len(ctypes)==1:
			BCS = CancerInfo.split('|')[3].split(':')[:4]
			BCS = [int(i) for i in BCS]
			BCS[i_ref] = 0 # setting REF to 0 as we only consider ALT
			MAX = max(BCS)
			index = np.argmax(BCS)
			BCS[index] = 0 # removing max to select next "max"
			MAX2 = max(BCS)

			ALT = 'ACTG'[index]
			BC = int(CancerInfo.split('|')[3].split(':')[index])
			CC = int(CancerInfo.split('|')[2].split(':')[index])
			VAF = round(BC/int(DP),4)
			MCF= round(CC/int(NC),4)

			FILTER = FILTER.replace('Multi-allelic,','')
			FILTER = FILTER.replace(',Multi-allelic','')
			FILTER = FILTER.replace('Multi-allelic','')

			if not(MAX2/MAX<0.05): # one alt 20x larger than the other)
				return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'Multi-Allelic'

			return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'PASS'
	else:	
		return ALT, FILTER, CTYPES, BC, CC, VAF, MCF,'PASS'
	
def BC_CC_filtering(STEP3FILTER,ALT,CancerInfo,min_ac_reads,min_ac_cells):
	alt_dict = {'A':0, 'C':1, 'T':2, 'G':3}
	i_alt = alt_dict[ALT[0]]
	try:
		CancerInfos = CancerInfo.split('|')
		BC = CancerInfos[3].split(':')[i_alt]
		CC = CancerInfos[2].split(':')[i_alt]
		if int(BC)<min_ac_reads or int(CC)<min_ac_cells:
			if STEP3FILTER == 'PASS':
				return 'LowDepth'
			else:
				return STEP3FILTER + ',LowDepth'
		else:
			return STEP3FILTER
	except AttributeError:
		if STEP3FILTER == 'PASS':
				return 'NoCov'
		else:
			return STEP3FILTER + ',NoCov'


def BetaBino_filtering(STEP3FILTER,CTYPES,CTYPESFILTER,POS):
	ctypes = CTYPES.split(',')
	if len(ctypes) == 1:
		if CTYPESFILTER in ['Non-Significant','Low-Significance']:
			if STEP3FILTER == 'PASS':
				return 'CancerNonSig'
			else:
				return STEP3FILTER + ',CancerNonSig'
		return STEP3FILTER
	elif len(ctypes) > 1:
		if ctypes[0] == 'Cancer':
			i_cancer = 0
			i_noncancer = 1
		elif ctypes[1] == 'Cancer':
			i_cancer = 1
			i_noncancer = 0
		if CTYPESFILTER.split(',')[i_cancer] in ['Non-Significant','Low-Significance']:
			if STEP3FILTER == 'PASS':
				return 'CancerNonSig'
			else:
				return STEP3FILTER + ',CancerNonSig'
		if CTYPESFILTER.split(',')[i_noncancer] in ['PASS','Low-Significance']:
			if STEP3FILTER == 'PASS':
				return 'NonCancerSig'
			else:
				return STEP3FILTER + ',NonCancerSig'
		return STEP3FILTER


def tag_clustered_SNVs(df, clust_dist):
	df2 = df[df['STEP3FILTER']=='PASS']
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
							   modify_filter(x['INDEX'],x['STEP3FILTER'], clust_dist, trash),
							   axis=1)
	print(updated_filter)
	return updated_filter

def modify_filter(INDEX, FILTER, clust_dist, trash):
	clustered = 'Clust_dist_{}'.format(str(clust_dist))
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


