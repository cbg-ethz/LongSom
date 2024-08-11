#!/usr/bin/python

import timeit
import argparse
import pandas as pd

def variant_calling_step3(file,outfile,filtered_outfile,deltaVAF,deltaCCF,cancer):

	#---
	# Command to focus only on high confidence cancer variants (HCCV)
	#---
	
	# Save comment lines
	f2 = open(outfile,'w')
	f3 = open(filtered_outfile,'w')
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
	f2.write('##INFO=HCCV_FILTER,Description=Filter status of the variant site for cell reannotation (high-confidence cancer variants)\n')		
	f3.write('##INFO=HCCV_FILTER,Description=Filter status of the variant site for cell reannotation (high-confidence cancer variants)\n')		
	f2.close()
	f3.close()
	
	input_df = pd.read_csv(file, sep='\t',comment='#',names=output_column_names)
	
	#Filtering homopolymeres, GnomAD/RNA editing/PON/clustered sites, 
	# multiallelic sites and sites with unsufficient # celltypes covered
	input_df = input_df[~input_df['FILTER'].str.contains('Min|LR|gnomAD|LC|RNA|llel', regex=True)]

	#Filtering short-read PON, except in chrM
	input_df = input_df[~((input_df['#CHROM']!='chrM') & (input_df['FILTER'].str.contains('SR')))]

	#Filtering loci passing all filters in cancer cells
	input_df = input_df[input_df['Cell_type_Filter'].str.contains('PASS', regex=True)]
	
	#Delta VAF and CCF filtering
	input_df['HCCV_FILTER'] = input_df.apply(lambda x: 
		HCCV_filtering(x['Cell_types'],x['Dp'],x['VAF'], x['CCF'],deltaVAF,deltaCCF,cancer), axis=1)

	filtered_df = input_df.copy(deep=True)[input_df['HCCV_FILTER'] == 'PASS']

	# Keeping only 1 position for clustered positions
	filtered_df['INDEX'] = filtered_df['#CHROM'] + ':' + filtered_df['Start'].astype(str) + ':' + filtered_df['ALT']
	idx = filtered_df['INDEX']
	a=[]
	for i in idx:
		chr,pos,base=i.split(':')
		a.append((chr,pos,base))

	b = sorted(a, key=lambda x: (x[0],x[1]))

	keep,trash=[],[]
	region = [0,0]
	for (chr,pos,base), (chr2,pos2,base2) in zip(b, b[1:]):
		if chr==chr2:
			if chr == 'chrM':
				keep.append(':'.join([chr,pos,base]))
				keep.append(':'.join([chr2,pos2,base2]))
			elif abs(int(pos)-int(pos2))<10000:
				print(region[1],int(pos))
				if region[1]==int(pos):
					region = [int(pos),int(pos2)]
					trash.append(':'.join([chr,pos,base]))
					trash.append(':'.join([chr2,pos2,base2]))
				else:
					keep.append(':'.join([chr,pos,base]))
					trash.append(':'.join([chr2,pos2,base2]))
					region = [int(pos),int(pos2)]
			else:
				keep.append(':'.join([chr,pos,base]))
		else:
			keep.append(':'.join([chr,pos,base]))
			keep.append(':'.join([chr2,pos2,base2]))
	keep = [ i for i in keep if i not in trash]
	filtered_df = filtered_df[filtered_df['INDEX'].isin(keep)]

	# Write outputs
	input_df.to_csv(outfile, sep='\t', index=False,  mode='a')
	filtered_df.to_csv(filtered_outfile, sep='\t', index=False,  mode='a')
	
def HCCV_filtering(CTYPES,DP,VAF,CCF,deltaVAFmin,deltaCCFmin, cancer):
	ctypes = CTYPES.split(',')
	
	#If only 1 cell type is called, no HCCV info
	if len(ctypes)==1:
		return 'NoInfo'
	else:
		DP1,DP2=DP.split(',')
		#Only look at higher depth loci
		if int(DP1)<10 or int(DP2)<10:
			return 'LowDepth'
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
			if VAF2>0.2:
				return 'HeterozygSite'
			# HCCF are variants with high VAF in cancer and low VAF in non-cancer cells
			elif deltaCCF < deltaCCFmin:
				return 'LowDeltaCCF'
			else:
				return 'PASS'
		#Same as above, but with inversed cell types.
		else:
			if VAF1>0.2:
				return 'HeterozygSite'
			elif -deltaCCF < deltaCCFmin:
				return 'LowDeltaCCF'
			else:
				return 'PASS'

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to perform the scRNA somatic variant calling')
	parser.add_argument('--infile', type=str, help='Input file with all samples merged in a single tsv', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--deltaVAF', type=float, default = 0.3, help='Delta VAF between cancer and non-cancer cells', required = True)
	parser.add_argument('--deltaCCF', type=float, default = 0.3, help='Delta CCF (cancer cell fraction) between cancer and non-cancer cells', required = True)
	parser.add_argument('--cancer_ctype', type=str, default = '', help='Name of cancer cell type in meta file', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	deltaVAF = args.deltaVAF
	deltaCCF = args.deltaCCF
	cancer_ctype = args.cancer_ctype

	# 1.2: Step 2: Add distance, editing and PoN filters
	print ('\n- High Confidence Cancer Variants calling\n')
	outfile3 = outfile + '.HCCV.unfiltered.tsv'
	filtered_outfile = outfile + '.HCCV.tsv'
	variant_calling_step3(infile,outfile3,filtered_outfile,deltaVAF,deltaCCF,cancer_ctype)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


