#!/usr/bin/python

import timeit
import argparse
import pandas as pd
import numpy as np
	
def variant_calling_step3(file,out_prefix,deltaVAF,deltaCCF,cancer,chrM_conta,min_ac_reads,min_ac_cells,clust_dist):

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
	input_df = input_df[input_df['Cell_types']!='NonCancer']

	# Filtering multiallelic sites, keeping only one allele if it's 50x larger than the other
	# Otherwise deleting the line (it messes with filters later on)
	# This is important as SComatic will often detect low expr secondary alt allele in very high expr. (e.g. in chrM)
	input_df[['ALT', 'FILTER', 'Cell_types', 'Bc', 'Cc', 'VAF', 'CCF','MultiAllelic_filter']] =  input_df.apply(lambda x: MultiAllelic_filtering(x['ALT'], x['FILTER'], x['Cell_types'],x['Bc'], x['Cc'], x['VAF'], x['CCF']), axis=1, result_type="expand") 
	input_df = input_df[input_df['MultiAllelic_filter']=='KEEP']
	input_df = input_df[output_column_names + ['INDEX']]

	#Special case for chrM due to contaminants:
	if chrM_conta == 'True':
		# Save chrM candidate SNVs to apply specific filters
		chrm_df = input_df[input_df['#CHROM']=='chrM'].copy()
		input_df = input_df[input_df['#CHROM']!='chrM']
		# Filtering homopolymeres, GnomAD/RNA editing/PON/clustered sites, 
		# multiallelic sites and sites with unsufficient # celltypes covered
		chrm_df = chrm_df[~chrm_df['FILTER'].str.contains('Min|LR|gnomAD|LC|RNA', regex=True)]
		if len(chrm_df)>0:
			# Applying deltaVAF and deltaCCF filters
			chrm_df['FILTER'] = chrm_df.apply(lambda x: 
				chrM_filtering(x['Cell_types'],x['Dp'],x['VAF'], x['CCF'],deltaVAF,deltaCCF,cancer), axis=1)
		else:
			chrM_conta = 'False'

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
	input_df = input_df[~input_df['Cell_type_Filter'].isin(['Non-Significant','Non-Significant'])]
	
	#Adding chrM SNVs detected above:
	if chrM_conta == 'True':
		unfiltered_df = pd.concat([input_df,chrm_df])
		unfiltered_df.to_csv(out_prefix+ '.calling.step3.unfiltered.tsv', sep='\t', index=False,  mode='a')
	else:
		input_df.to_csv(out_prefix+ '.calling.step3.unfiltered.tsv', sep='\t', index=False,  mode='a')

	# Filter 4: Noise filter:
	input_df = input_df[~input_df['FILTER'].str.contains('Noisy_site')] #multi allelic sites are filtered above
	
	# Filter 5: Homopolymer filter:
	input_df = input_df[~input_df['FILTER'].str.contains('LC_Upstream|LC_Downstream', regex = True)] #multi allelic sites are filtered above
	
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
	if chrM_conta == 'True':
		input_df = pd.concat([input_df,chrm_df])

	#Filtering PASS SNVs (only left in chrM)
	filtered_df = input_df[input_df['FILTER'] =='PASS']

	# Write output
	filtered_df.to_csv(out_prefix+ '.calling.step3.tsv', sep='\t', index=False,  mode='a')
	
def chrM_filtering(CTYPES,DP,VAF,CCF,deltaVAFmin,deltaCCFmin,cancer):
	ctypes = CTYPES.split(',')
	VAF=VAF
	CCF=CCF
	DP=DP
	if len(ctypes)>1:
		DP1,DP2=DP.split(',')
		#Only look at higher depth loci
		if int(DP1)<100 or int(DP2)<100: 
			return 'LowDepth'
		else: 
			VAFs = VAF.split(',')
			VAF1 = float(VAFs[0])
			VAF2 = float(VAFs[1])
			deltaVAF = VAF1-VAF2
			CCFs = CCF.split(',')
			CCF1 = float(CCFs[0])
			CCF2 = float(CCFs[1])
			deltaCCF = CCF1-CCF2
			if ctypes[0]==cancer:
				#If Non-cancer VAF is high, it's likely an heterozygous site
				if deltaVAF < deltaVAFmin:
						return 'LowDeltaVAF'
				# HCCF are variants with high VAF in cancer and low VAF in non-cancer cells
				elif deltaCCF < deltaCCFmin:
						return 'LowDeltaCCF'
				else:
						return 'PASS'
			#Same as above, but with inversed cell types.
			else:
				if -deltaVAF < deltaVAFmin:
						return 'LowDeltaVAF'
				elif -deltaCCF < deltaCCFmin:
						return 'LowDeltaCCF'
				else:
						return 'PASS'
	else:
		if ctypes[0]!=cancer:
			return 'NonCancer'
		elif int(DP)<50:
			return 'LowDepth'
		elif float(VAF)<0.1:
			return 'LowVAF'
		elif float(CCF)<0.1:
			return 'LowCCF'
		else:
			return 'PASS'

def MultiAllelic_filtering(ALT, FILTER, CTYPES, BC, CC, VAF, CCF):

	if len(ALT.split('|'))>1:
		if len(CTYPES.split(','))>1:
			if len(BC.split(',')[1].split('|'))>1:
				BCS = [int(i) for i in BC.split(',')[1].split('|')]
				MAX = max(BCS)
				index = np.argmax(BCS)
				BCS[index] = 0 # removing max to select next "max"
				MAX2 = max(BCS)
				if not(MAX2/MAX<0.02): # one alt 50x larger than the other)
					return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
			else:
				return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
				
			ALT_Cancer = ALT.split(',')[1].split('|')[index]
			BC_Cancer = BC.split(',')[1].split('|')[index]
			CC_Cancer = CC.split(',')[1].split('|')[index]
			VAF_Cancer = VAF.split(',')[1].split('|')[index]
			CCF_Cancer = CCF.split(',')[1].split('|')[index]
			try:
				ALT_NonCancer = ALT.split(',')[0].split('|')[index]
				BC_NonCancer = BC.split(',')[0].split('|')[index]
				CC_NonCancer = CC.split(',')[0].split('|')[index]
				VAF_NonCancer = VAF.split(',')[0].split('|')[index]
				CCF_NonCancer = CCF.split(',')[0].split('|')[index]
			except:
				ALT_NonCancer = ALT.split(',')[0].split('|')[0]
				BC_NonCancer = BC.split(',')[0].split('|')[0]
				CC_NonCancer = CC.split(',')[0].split('|')[0]
				VAF_NonCancer = VAF.split(',')[0].split('|')[0]
				CCF_NonCancer = CCF.split(',')[0].split('|')[0]
			ALT = ','.join([ALT_NonCancer,ALT_Cancer])
			BC = ','.join([BC_NonCancer,BC_Cancer])
			CC = ','.join([CC_NonCancer,CC_Cancer])
			VAF = ','.join([VAF_NonCancer,VAF_Cancer])
			CCF = ','.join([CCF_NonCancer,CCF_Cancer])
			return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'KEEP'
		else:
			BCS = [int(i) for i in BC.split('|')]
			if min(BCS)/max(BCS)<0.02: # one alt 50x larger than the other)
				index = np.argmax(BCS)
			else:
				return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
				
			ALT = ALT.split('|')[index]
			BC = BC.split('|')[index]
			CC = CC.split('|')[index]
			VAF = VAF.split('|')[index]
			CCF = CCF.split('|')[index]
			return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'KEEP'
	else:	
		if 'Multi-allelic' in FILTER: #Dissonant alleles
			return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
		return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'KEEP'
	
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
	parser.add_argument('--deltaCCF', type=float, default = 0.3, help='Delta CCF (cancer cell fraction) between cancer and non-cancer cells', required = True)
	parser.add_argument('--cancer_ctype', type=str, default = '', help='Name of cancer cell type in meta file', required = False)
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
	deltaCCF = args.deltaCCF
	cancer_ctype = args.cancer_ctype
	chrM_conta = args.chrM_contaminant
	min_ac_reads = args.min_ac_reads
	min_ac_cells = args.min_ac_cells
	clust_dist = args.clust_dist

	# 1.1: Step 3: Filter only PASS mutations, not within clust_dist of each other, and apply specific chrM filters
	print ('\n- Variant calling step 3\n')
	variant_calling_step3(infile,out_prefix,deltaVAF,deltaCCF,cancer_ctype,chrM_conta,min_ac_reads,min_ac_cells,clust_dist)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


