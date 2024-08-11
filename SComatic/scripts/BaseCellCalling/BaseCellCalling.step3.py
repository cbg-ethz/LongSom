#!/usr/bin/python

import timeit
import argparse
import pandas as pd
	
def variant_calling_step3(file,out_prefix,deltaVAF,deltaCCF,cancer,chrM_conta,min_ac_reads,clust_dist,PoN_LR):

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
	input_df['INDEX'] = input_df['#CHROM'] + ':' + input_df['Start'].astype(str) + ':' + input_df['ALT'].str.split(',', n=1, expand=True)[0]
	#Filtering multiallelic sites
	input_df = input_df[~input_df['VAF'].str.contains('\|')]	

	#Special case for chrM due to contaminants:
	if chrM_conta == 'True':
		# Save chrM candidate SNVs to apply specific filters
		chrm_df = input_df[input_df['#CHROM']=='chrM'].copy()
		input_df = input_df[input_df['#CHROM']!='chrM']
		# Filtering homopolymeres, GnomAD/RNA editing/PON/clustered sites, 
		# multiallelic sites and sites with unsufficient # celltypes covered
		chrm_df = chrm_df[~chrm_df['FILTER'].str.contains('Min|LR|gnomAD|LC|RNA|llel', regex=True)]
		# Applying deltaVAF and deltaCCF filters
		chrm_df['FINAL_FILTER'] = chrm_df.apply(lambda x: 
			chrM_filtering(x['Cell_types'],x['Dp'],x['VAF'], x['CCF'],deltaVAF,deltaCCF,cancer), axis=1)
		#Only keep one ALT base in case it's duplicated 
		input_df['ALT'] = input_df['ALT'].str.split(',').str[0]
	
	#Filter mutations found in cancer cells
	input_df = input_df[input_df['Cell_types'] == cancer]
	#Filtering min alt reads:
	input_df = input_df[input_df['Cc'].astype(int)>=min_ac_reads]


	#Filtering out SNVs close to each other (likely mismapping or CNV event)
	input_df['FINAL_FILTER'] = tag_clustered_SNVs(input_df, clust_dist)
	
	#Special case for chrM due to contaminants:
	if chrM_conta == 'True':
		input_df = pd.concat([input_df,chrm_df])

	input_df.to_csv(out_prefix+ '.calling.step3.unfiltered.tsv', sep='\t', index=False,  mode='a')

	#Filtering PASS SNVs
	if PoN_LR!='True':
		filtered_df = input_df[input_df['FINAL_FILTER'].isin(['PASS','Clustered','PoN_LR','Clustered,PoN_LR'])]
	else:
		filtered_df = input_df[input_df['FINAL_FILTER'].isin(['PASS','Clustered'])]

	# Write output
	filtered_df.to_csv(out_prefix+ '.calling.step3.tsv', sep='\t', index=False,  mode='a')
	
def chrM_filtering(CTYPES,DP,VAF,CCF,deltaVAFmin,deltaCCFmin,cancer):
	ctypes = CTYPES.split(',')
	VAF=VAF.split('|')[0]
	CCF=CCF.split('|')[0]
	DP=DP.split('|')[0]
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
		else:
			return 'PASS'
	
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

	df['FINAL_FILTER'] = df.apply(lambda x : modify_filter(x['INDEX'],x['FILTER'], clust_dist, trash))
	return df['FINAL_FILTER']

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
	parser.add_argument('--PoN_LR', type=str, default = 'True', help='Filter mutations based on PoN_LR (default: True)', required = False)
	parser.add_argument('--min_ac_reads', type=int, default = 5, help='Minimum ALT reads', required = False)
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
	PoN_LR = args.PoN_LR
	min_ac_reads = args.min_ac_reads
	clust_dist = args.clust_dist

	# 1. Variant calling
	print ('\n------------------------------')
	print ('Variant calling')
	print ('------------------------------\n')

	# 1.2: Step 2: Add distance, editing and PoN filters
	print ('\n- Variant calling step 3\n')
	variant_calling_step3(infile,out_prefix,deltaVAF,deltaCCF,cancer_ctype,chrM_conta,min_ac_reads,clust_dist,PoN_LR)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


