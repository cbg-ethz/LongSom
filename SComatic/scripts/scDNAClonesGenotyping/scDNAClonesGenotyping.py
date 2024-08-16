#!/usr/bin/python

import timeit
import argparse
import pandas as pd
import pysam
import timeit
import multiprocessing as mp
import argparse
import glob
import os
import math

def collect_result(results):
	pass

def BaseCount(LIST, REF_BASE):
	Bases=['A','C','T','G','N']
	
	# Dictinary with base counts
	NUCLEOTIDES = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'I': 0, 'D' : 0, 'N': 0, 'R' : 0, 'Other' : 0}
	
	# Update dictionary counts
	for x in LIST:
		if x.upper() in Bases:
			NUCLEOTIDES[x.upper()] += 1
		elif '-' in x:
			NUCLEOTIDES['D'] += 1
		elif '+' in x:
			NUCLEOTIDES['I'] += 1
		else:
			NUCLEOTIDES['Other'] += 1
	
	# Calculate Alternative count of the most frequent alternative allele	
	Alternative = ['A','C','T','G','I', 'D']
	
	ALTERNATIVE = {k:v for k,v in NUCLEOTIDES.items() if k != REF_BASE.upper() and k in Alternative}

	AC = max(ALTERNATIVE.items(), key=lambda x: x[1])[1]
	
	if (AC > 0):
		listOfKeys = list()
		# Iterate over all the items in dictionary to find keys with max value
		for key, value in ALTERNATIVE.items():
			if value == AC:
				listOfKeys.append(key)
		
		MAX_ALT = ','.join(sorted(listOfKeys))
	else:
		MAX_ALT = '.'
		
	return ([NUCLEOTIDES,MAX_ALT,AC])

#@profile
def EasyReadPileup(LIST, REF_BASE):
	Bases = set(['A','C','T','G','N'])
	
	AC = 0

	# Create a new list based on pileup reading
	NEW_LIST = []
	for x in LIST:
		LEN = len(x)
		UPPER = x.upper()
		if UPPER in Bases:
			NEW_LIST.append(UPPER)
			# Check for Alt counts
			if (UPPER != REF_BASE):
				AC = AC + 1
		elif LEN > 1 and x[1] == '-':
			D = 'D'
			NEW_LIST.append(D)
			AC = AC + 1
		elif LEN > 1 and x[1] == '+':
			I = 'I'
			NEW_LIST.append(I)
			AC = AC + 1
		elif x == '*':
			NEW_LIST.append('O')
		else:
			NEW_LIST.append('NA')
			
	return (NEW_LIST,AC)

def run_interval(code,var_dict,meta_dict,clones,majorclones,BAM,FASTA,tmp_dir,BQ,MQ,ALT_FLAG):
	# List of positions	
	interval = var_dict[code]
	
	# Create a dictionary with the positions to analyse
	CHROM = []
	Target_sites = {}
	for group in interval:
		group = group.split("\t")

		# Chrom
		# Extract as well info from windows
		CHROM.append(group[0])

		# Position	
		# Consider that pysam uses 0-based coordinates (REAL_POS - 1)
		pos_temp = int(group[1])
		pos_temp = pos_temp - 1

		# Fill sub-dictionary with variant information
		REF = group[3]
		ALT = group[4].split(',')[0]
		CELL_TYPE_EXPECTED = group[6]
		NUM_CELLS_WITH_MUTATION_EXPECTED = group[13]

		Target_sites[pos_temp] = [REF,ALT,CELL_TYPE_EXPECTED,NUM_CELLS_WITH_MUTATION_EXPECTED]

	# Coordinates to analyse
	List_of_target_sites = set(Target_sites.keys())
	CHROM  =  CHROM[0]
	START = min(List_of_target_sites) - 1
	END = max(List_of_target_sites) + 1
	
	# Temporary out file
	ID = '_'.join([str(CHROM), str(min(List_of_target_sites)), str(max(List_of_target_sites))])
	out_temp = tmp_dir + '/' + ID + '.CloneGenotype.temp'
	outfile = open(out_temp,'w')
	# Get pileup read counts from coordinates
	bam = pysam.AlignmentFile(BAM)
	i = bam.pileup(CHROM, START, END, min_base_quality = BQ, min_mapping_quality = MQ, ignore_overlaps = False, max_depth = 200000)
	
	# Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
	inFasta = pysam.FastaFile(FASTA)
	
	#Cell were the barcode is found
	CLONES = {pos:{clone:{'Dp':0,'Alt':0} for clone in clones} for pos in List_of_target_sites }
	MAJORCLONES = {pos:{clone:{'Dp':0,'Alt':0} for clone in majorclones} for pos in List_of_target_sites }
	
	# Run it for each position in pileup
	for p in i:
		POS=p.pos
		if POS in List_of_target_sites:
			# Expected values. Variants to check
			Ref_exp,Alt_exp,Cell_type_exp,Num_cells_exp = Target_sites[POS]

			# Get reference base from fasta file
			ref_base = inFasta.fetch(CHROM, POS, POS+1)
			
			# Get pileup info
			PILEUP_LIST = p.get_query_sequences(mark_matches=True, add_indels=True)
			NEW_PILEUP_LIST,AC = EasyReadPileup(PILEUP_LIST, ref_base)

			# Expected bases
			if (ALT_FLAG == 'All'):
				Bases = set(['A','C','T','G','I','D','N'])
				idx_reads = [x for x in range(0,len(NEW_PILEUP_LIST)) if NEW_PILEUP_LIST[x] in Bases]
			else:
				idx_reads = [x for x in range(0,len(NEW_PILEUP_LIST)) if NEW_PILEUP_LIST[x] == Alt_exp]

			# Get reads info
			reads = p.pileups

			for read_i in idx_reads:
				
				read_alignment = reads[read_i].alignment

				# Alt observed in the read
				Alt_obs = NEW_PILEUP_LIST[read_i]
				try:
					barcode = read_alignment.opt("CB")
					clone = meta_dict[barcode]
					majclone = clone.rsplit('_',1)[0]
				except:
					continue

				# Remove low quality reads: Seconday alignments, duplicate reads, supplementary alignments
				if read_alignment.is_secondary == False and read_alignment.is_duplicate == False and read_alignment.is_supplementary == False:
					#compute depth and #alt reads per cell
					if str(Alt_obs) == str(Alt_exp):
						CLONES[POS][clone]['Dp']+=1
						CLONES[POS][clone]['Alt']+=1
						MAJORCLONES[POS][majclone]['Dp']+=1
						MAJORCLONES[POS][majclone]['Alt']+=1
					else:
						CLONES[POS][clone]['Dp']+=1
						MAJORCLONES[POS][majclone]['Dp']+=1

	#Clones and major clones VAF
	for POS in CLONES.keys():
		for CLONE in CLONES[POS].keys():
			#Initialize values
			DP = CLONES[POS][CLONE]['Dp']
			ALT = CLONES[POS][CLONE]['Alt'] 
			CLONES[POS][CLONE]['VAF']  = '.'
			if DP > 0:
				CLONES[POS][CLONE]['VAF'] = round(ALT/DP,4)
		for CLONE in MAJORCLONES[POS].keys():
			#Initialize values
			DP = MAJORCLONES[POS][CLONE]['Dp']
			ALT = MAJORCLONES[POS][CLONE]['Alt'] 
			MAJORCLONES[POS][CLONE]['VAF'] = '.'
			if DP > 0:
				MAJORCLONES[POS][CLONE]['VAF'] = round(ALT/DP,4)

		Ref_exp,Alt_exp,Cell_type_exp,Num_cells_exp = Target_sites[POS]
		INDEX = str(CHROM) + ':' + str(POS+1) + ':' + Alt_exp.split(',')[0]
		Group = [str(CHROM), str(POS+1),str(POS+1),Ref_exp,Alt_exp,INDEX]  #CTYPE,str(DP),str(ALT),str(VAF),str(BETABIN),str(MUTATED),str(BINARIZED),INDEX]
		for CLONE in MAJORCLONES[POS]:
			Group.append(str(MAJORCLONES[POS][CLONE]['Dp']))
			Group.append(str(MAJORCLONES[POS][CLONE]['Alt']))
			Group.append(str(MAJORCLONES[POS][CLONE]['VAF']))
		for CLONE in CLONES[POS]:
			Group.append(str(CLONES[POS][CLONE]['Dp']))
			Group.append(str(CLONES[POS][CLONE]['Alt']))
			Group.append(str(CLONES[POS][CLONE]['VAF']))
		Group = '\t'.join(Group) + '\n'
		outfile.write(Group)		
			
	inFasta.close()
	bam.close()
	outfile.close()


def meta_to_dict(meta_file,tissue,ctype):
	metadata = pd.read_csv(meta_file, delimiter = "\t")

	# Clean index column
	#metadata['Index_clean'] = metadata['Index'].str.replace('-.*$','',regex=True)

	# If tissue provided, append tissue id to cell type to be recognised in downstream analysis
	if tissue == None:
		metadata['Cell_type_clean'] = metadata[ctype].str.replace(' ','_',regex=True)
	else:
		tissue = tissue.replace(" ", "_")
		metadata['Cell_type_clean'] = metadata[ctype].str.replace(' ','_',regex=True)
		metadata['Cell_type_clean'] = str(tissue) + '__' + metadata['Cell_type_clean'].astype(str)

	# Create dicitionary with cell types and cell barcodes
	DICT = dict(zip( metadata['Index'], metadata['Cell_type_clean']))
	CLONES = list(metadata['Cell_type_clean'].unique())
	MAJORCLONES = list(set([c.rsplit('_',1)[0] for c in CLONES]))

	del metadata
	return(DICT, CLONES, MAJORCLONES)


def build_dict_variants(variant_file,window):
	DICT_variants = {}
	with open(variant_file, 'r') as f:
		for line in f:
			if not line.startswith('#') and not line.startswith('Chr'):
				line = line.rstrip('\n')
				elements = line.split('\t')

				# Chromosome
				CHROM = elements[0]

				# Position
				POS = int(elements[1])
				WIND = math.floor(POS / float(window))

				CODE = CHROM + '_' + str(WIND)
				if not CODE in DICT_variants.keys():
					DICT_variants[CODE] = [line]
				else:
					DICT_variants[CODE].append(line)

	return (DICT_variants)

def concatenate_sort_clonetemp_files_and_write(out_prefix, tmp_dir, clones, majorclones,):
	# Get the file paths
	all_files = glob.glob(tmp_dir + '/*.CloneGenotype.temp')
	clones = majorclones + clones
	# Load as panda files
	if (len(all_files) > 0):
		
		## Organize files in dictionaries of chromosomes and starts
		Dictionary_of_files = {}
		for filename in all_files:
			basename = os.path.basename(filename)
			
			# It can consider not standard chromosome nomenglatures (with dots)
			coordinates = basename.replace('.CloneGenotype.temp','')

			CHROM, START, END = coordinates.rsplit("_",2)
			
			START = int(START)
			
			if (CHROM not in Dictionary_of_files):
				Dictionary_of_files[CHROM] = {}
				Dictionary_of_files[CHROM][START] = filename
			else:
				Dictionary_of_files[CHROM][START] = filename
		

		## Write in the final output file
		out = open(out_prefix+'.CloneGenotype.tsv','w')
		Header=['#CHROM','Start', 'End', 'REF', 'ALT_expected', 'INDEX']
		for clone in clones:
			Header.append(clone + '_DP')
			Header.append(clone + '_ALT')
			Header.append(clone + '_VAF')
		out.write('\t'.join(Header)+'\n')
		
		# Move thrugh filenanes by sorted coordinates
		for chrom in sorted(Dictionary_of_files.keys()):
			for start in sorted(Dictionary_of_files[chrom].keys()):
				filename = Dictionary_of_files[chrom][start]
				
				with open(filename, 'r') as f:
					out.write(f.read())
		
				# Remove temp file
				os.remove(filename)
		out.close()
			
	
	else:
		# If temp files not found, print message
		print ('No temporary files found')

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get the alleles observed in each unique cell for the variant sites')
	parser.add_argument('--bam', type=str, default=1, help='Tumor bam file to be analysed', required = True)
	parser.add_argument('--infile', type=str, default=1, help='Base calling file (obtained by BaseCellCalling.step2.py), ideally only the PASS variants', required = True)
	parser.add_argument('--ref', type=str, default=1, help='Reference genome. *fai must be available in the same folder as reference', required = True)
	parser.add_argument('--meta', type=str, default=1, help='Metadata with cell barcodes per cell type', required = True)
	parser.add_argument('--fusions', type=str, help='Fusions file from CTAT_fusion', required = False)
	parser.add_argument('--outfile', default = 'Matrix.tsv', help='Out file', required = False)
	parser.add_argument('--alt_flag', default = 'All',choices = ['Alt','All'], help='Flag to search for cells carrying the expected alt variant (Alt) or all cells independent of the alt allele observed (All)', required = False)
	parser.add_argument('--nprocs',default = 1, help='Number of processes [Default = 1]',required=False,type = int)
	parser.add_argument('--bin', type=int, default=50000, help='Bin size for running the analysis [Default 50000]', required = False)
	parser.add_argument('--min_bq', type=int, default = 30, help='Minimum base quality permited for the base counts. Default = 30', required = False)
	parser.add_argument('--min_mq', type=int, default = 255, help='Minimum mapping quality required to analyse read. Default = 255', required = False)
	parser.add_argument('--tissue', type=str, default=None, help='Tissue of the sample', required = False)
	parser.add_argument('--ctype', type=str, default='Cell_type', help='Name of cell type column in meta file', required = False)
	parser.add_argument('--tmp_dir', type=str, default = 'tmpDir', help='Temporary folder for tmp files', required = False)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	BAM = args.bam
	FASTA = args.ref
	CORE = args.nprocs
	out_prefix = args.outfile
	window = args.bin
	tmp_dir = args.tmp_dir
	meta_file = args.meta
	variant_file = args.infile
	tissue = args.tissue
	BQ = args.min_bq
	MQ = args.min_mq
	ALT_FLAG = args.alt_flag
	ctype = args.ctype

	# Set outfile name
	print("Outfile prefix: " , out_prefix ,  "\n") 

	# 1. Check if temp dir is created
	if (tmp_dir != '.'):
		try:
			# Create target Directory
			os.mkdir(tmp_dir)
			print("Directory " , tmp_dir ,  " created\n") 
		except FileExistsError:
			print("Directory " , tmp_dir ,  " already exists\n")
	else:
		print("Not temp directory specified, using working directory as temp") 

	# 2. Create dictionaries with variants and cell barcodes
	META_DICT, CLONES, MAJORCLONES = meta_to_dict(meta_file,tissue,ctype)
	DICT_VARIANTS = build_dict_variants(variant_file,window)
	
	# 3. Code to run in parallel all bins
	if (CORE > 1):
		pool = mp.Pool(CORE)
		
		# Step 3.1: Use loop to parallelize
		for row in DICT_VARIANTS.keys():
			# This funtion writes in temp files the results
			pool.apply_async(run_interval, args=(row,DICT_VARIANTS,META_DICT,CLONES,MAJORCLONES,BAM,FASTA,tmp_dir,BQ,MQ,ALT_FLAG), callback=collect_result)
				   
		# Step 3.2: Close Pool and let all the processes complete	
		pool.close()
		pool.join()
	else:
		for row in DICT_VARIANTS.keys():
			# This funtion writes in temp files the results
			collect_result(run_interval(row,DICT_VARIANTS,META_DICT,CLONES,MAJORCLONES,BAM,FASTA,tmp_dir,BQ,MQ,ALT_FLAG))
	
	# 4. Write "Long" file
	concatenate_sort_clonetemp_files_and_write(out_prefix, tmp_dir, CLONES, MAJORCLONES)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')

