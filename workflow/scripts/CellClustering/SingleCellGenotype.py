import pysam
import timeit
import multiprocessing as mp
import argparse
import pandas as pd
import glob
import os
import math
import numpy as np
from natsort import natsorted
from scipy.stats import betabinom

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

def run_interval(code,var_dict,meta_dict,BAM,FASTA,tmp_dir,BQ,MQ,ALT_FLAG,alpha2,beta2,pval,chrm_conta):
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
	out_temp = tmp_dir + '/' + ID + '.SingleCellCounts.temp'
	outfile = open(out_temp,'w')

	# Get pileup read counts from coordinates
	bam = pysam.AlignmentFile(BAM)
	i = bam.pileup(CHROM, START, END, min_base_quality = BQ, min_mapping_quality = MQ, ignore_overlaps = False, max_depth = 200000)
	
	# Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
	inFasta = pysam.FastaFile(FASTA)
	
	#Cell were the barcode is found
	CELLS = {pos:{barcode:{'Dp':0,'Alt':0} for barcode in meta_dict.keys()} for pos in List_of_target_sites }
	
	# Run it for each position in pileup
	for p in i:
		POS=p.pos
		if POS in List_of_target_sites:
			# Expected values. Variants to check
			Ref_exp,Alt_exp,Cell_type_exp,Num_cells_exp = Target_sites[POS]

			# Get reference base from fasta file
			ref_base = inFasta.fetch(CHROM, POS, POS+1)
			
			# Get pileup info
			READS = p.get_query_names()
			QUALITIES = p. get_query_qualities()
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
					barcode = barcode.split("-")[0]
					ctype = meta_dict[barcode]
				except:
					continue

				# Remove low quality reads: Seconday alignments, duplicate reads, supplementary alignments
				if read_alignment.is_secondary == False and read_alignment.is_duplicate == False and read_alignment.is_supplementary == False:
					#compute depth and #alt reads per cell
					if str(Alt_obs) == str(Alt_exp):
						CELLS[POS][barcode]['Dp']+=1
						CELLS[POS][barcode]['Alt']+=1
					else:
						CELLS[POS][barcode]['Dp']+=1

	#Compute VAF and establish if the locus is mutated in each cell with coverage
	for POS in CELLS.keys():
		for bc in CELLS[POS].keys():
			try:
				CTYPE = meta_dict[bc]
			except:
				continue
			#Initialize values
			DP = CELLS[POS][bc]['Dp']
			ALT = CELLS[POS][bc]['Alt'] 
			VAF = '.'
			BETABIN = '.'
			MUTATED = 'NoCoverage'
			if DP > 0:
				VAF = round(ALT/DP,4)
				if ALT > 0:
					#Special case for chrM due to contaminants:
					if chrm_conta == 'True' and str(CHROM) == 'chrM':
							if VAF < 0.3:
								MUTATED = 'LowVAFChrM'
							else:
								MUTATED = 'PASS'
					#All non-chrM chromosomes
					else:
						BETABIN = round(betabinom.sf(ALT-0.001, DP, alpha2, beta2),4)
						#test if the betabin p-value is significant
						if BETABIN < pval:
							MUTATED = 'PASS'
						else:
							MUTATED = 'BetaBin_problem'	
				else:
					MUTATED = 'NoAltReads'

			if MUTATED == "PASS":
				BINARIZED = 1
			elif MUTATED == "NoCoverage":
				BINARIZED = 3
			else:
				BINARIZED = 0

			Ref_exp,Alt_exp,Cell_type_exp,Num_cells_exp = Target_sites[POS]
			INDEX = str(CHROM) + ':' + str(POS+1) + ':' + Alt_exp.split(',')[0]
			Group = [str(CHROM), str(POS+1),str(POS+1),Ref_exp,Alt_exp,str(Cell_type_exp),str(Num_cells_exp),bc,CTYPE,str(DP),str(ALT),str(VAF),str(BETABIN),str(MUTATED),str(BINARIZED),INDEX]
			Group = '\t'.join(Group) + '\n'
			outfile.write(Group)
			
	inFasta.close()
	bam.close()
	outfile.close()

def meta_to_dict(meta_file,tissue):
	metadata = pd.read_csv(meta_file, delimiter = "\t")

	# Clean index column
	metadata['Index_clean'] = metadata['Index'].str.replace('-.*$','',regex=True)

	# If tissue provided, append tissue id to cell type to be recognised in downstream analysis
	if tissue == None:
		metadata['Cell_type_clean'] = metadata['Cell_type'].str.replace(' ','_',regex=True)
	else:
		tissue = tissue.replace(" ", "_")
		metadata['Cell_type_clean'] = metadata['Cell_type'].str.replace(' ','_',regex=True)
		metadata['Cell_type_clean'] = str(tissue) + '__' + metadata['Cell_type_clean'].astype(str)

	# Create dicitionary with cell types and cell barcodes
	DICT = metadata.set_index('Index_clean')['Cell_type_clean'].to_dict()
	ALL_CELL_TYPES = metadata['Cell_type_clean'].unique()

	del metadata

	return(DICT, ALL_CELL_TYPES)


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


def concatenate_sort_temp_files_and_write(out_prefix, tmp_dir):
	# Get the file paths
	all_files = glob.glob(tmp_dir + '/*.SingleCellCounts.temp')
	
	# Load as panda files
	if (len(all_files) > 0):
		
		## Organize files in dictionaries of chromosomes and starts
		Dictionary_of_files = {}
		for filename in all_files:
			basename = os.path.basename(filename)
			
			# It can consider not standard chromosome nomenglatures (with dots)
			coordinates = basename.replace('.SingleCellCounts.temp','')

			CHROM, START, END = coordinates.rsplit("_",2)
			
			START = int(START)
			
			if (CHROM not in Dictionary_of_files):
				Dictionary_of_files[CHROM] = {}
				Dictionary_of_files[CHROM][START] = filename
			else:
				Dictionary_of_files[CHROM][START] = filename
		

		## Write in the final output file
		out = open(out_prefix+'.SingleCellGenotype.tsv','w')
		Header=['#CHROM','Start', 'End', 'REF', 'ALT_expected', 'Cell_type_expected','Num_cells_expected','CB','Cell_type_observed','Dp','ALT','VAF','BetaBin','MutationStatus','BinMutationStatus','INDEX']
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

def collect_cells_with_fusions(fusion_file):
	fusions = pd.read_csv(fusion_file, sep='\t')

	# Create Fusion:barcode index to remove duplicates
	fusions['INDEX'] = fusions['#FusionName'] + ':' + fusions['BC']

	# Drop duplicates
	fusions = fusions.drop_duplicates(subset='INDEX', keep="last")

	# Compute fusions to append them to SingleCellGenotype.tsv
	lines = []
	for i,row in fusions.iterrows():
		l = ['.','.','.','.','.','.','.',row['BC'],'.', 1,1,1,'.','.',1,'zzz:'+row['#FusionName']]
		lines.append(l)

	return pd.DataFrame(lines)

def sort_chr_index(df):
	renamed = [i.replace('chrM','chrZ') for i in df.index]
	df.index = renamed
	df = df.reindex(natsorted(df.index))
	renamed = [i.replace('chrZ','chrM').replace('zzz:','') for i in df.index]
	df.index = renamed 
	return df


def pivot_long_dataframe(out_prefix, fusions):
	#Pivot long SingleCellGenotype file into wide cell-SNV matrixes
	long_df = pd.read_csv(out_prefix + '.SingleCellGenotype.tsv', sep = '\t')
	
	if not fusions.empty:
	#Concatenate fusions with SingleCellGenotype.tsv
		fusions.columns = long_df.columns
		long_df = pd.concat([long_df,fusions]).fillna(3)

	# Cell-SNV matrix with reads depth as values
	Dp_df = long_df.pivot(index='INDEX', columns='CB', values='Dp')
	Dp_df = sort_chr_index(Dp_df)
	Dp_df.to_csv(out_prefix + '.DpMatrix.tsv', sep='\t', index=True)

	# Cell-SNV matrix with raw ALT reads counts as values
	Alt_df = long_df.pivot(index='INDEX', columns='CB', values='ALT')
	Alt_df = sort_chr_index(Alt_df)
	Alt_df.to_csv(out_prefix + '.AltMatrix.tsv', sep='\t', index=True)

	# Cell-SNV matrix with VAFs as values
	VAF_df = long_df.pivot(index='INDEX', columns='CB', values='VAF')
	VAF_df = sort_chr_index(VAF_df)
	VAF_df.to_csv(out_prefix + '.VAFMatrix.tsv', sep='\t', index=True)

	# Cell-SNV matrix with binarized mutation status 
	# (1: mutated, 0: non-mutated, 3: no coverage)
	Bin_df = long_df.pivot(index='INDEX', columns='CB', values='BinMutationStatus')
	Bin_df = sort_chr_index(Bin_df)
	Bin_df.to_csv(out_prefix + '.BinaryMatrix.tsv', sep='\t', index=True)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get the SNV/fusions observed in each unique cell')
	parser.add_argument('--bam', type=str, default=1, help='Tumor bam file to be analysed', required = True)
	parser.add_argument('--infile', type=str, default=1, help='Base calling file (obtained by BaseCellCalling.step2.py), ideally only the PASS variants', required = True)
	parser.add_argument('--ref', type=str, default=1, help='Reference genome. *fai must be available in the same folder as reference', required = True)
	parser.add_argument('--meta', type=str, default=1, help='Metadata with cell barcodes per cell type', required = True)
	parser.add_argument('--fusions', type=str, help='Fusions file from CTAT_fusion', nargs='?', const='', required = True)
	parser.add_argument('--outfile', default = 'Matrix.tsv', help='Out file', required = False)
	parser.add_argument('--alt_flag', default = 'All',choices = ['Alt','All'], help='Flag to search for cells carrying the expected alt variant (Alt) or all cells independent of the alt allele observed (All)', required = False)
	parser.add_argument('--nprocs',default = 1, help='Number of processes [Default = 1]',required=False,type = int)
	parser.add_argument('--bin', type=int, default=50000, help='Bin size for running the analysis [Default 50000]', required = False)
	parser.add_argument('--min_bq', type=int, default = 30, help='Minimum base quality permited for the base counts. Default = 30', required = False)
	parser.add_argument('--min_mq', type=int, default = 255, help='Minimum mapping quality required to analyse read. Default = 255', required = False)
	parser.add_argument('--tissue', type=str, default=None, help='Tissue of the sample', required = False)
	parser.add_argument('--tmp_dir', type=str, default = 'tmpDir', help='Temporary folder for tmp files', required = False)
	parser.add_argument('--alpha2', type=float, default = 0.2474528917555431, help='Alpha parameter for Beta-binomial distribution of read counts. [Default: 0.260288007167716]', required = False)
	parser.add_argument('--beta2', type=float, default = 162.03696139428595, help='Beta parameter for Beta-binomial distribution of read counts. [Default: 173.94711910763732]', required = False)
	parser.add_argument('--pvalue', type=float, default = 0.01, help='P-value for the beta-binomial test to be significant', required = False)
	parser.add_argument('--chrM_contaminant', type=str, default = 'True', help='Use this option if chrM contaminants are observed in non-cancer cells', required = False)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	BAM = args.bam
	FASTA = args.ref
	CORE = args.nprocs
	fusion_file = args.fusions
	out_prefix = args.outfile
	window = args.bin
	tmp_dir = args.tmp_dir
	meta_file = args.meta
	variant_file = args.infile
	tissue = args.tissue
	BQ = args.min_bq
	MQ = args.min_mq
	ALT_FLAG = args.alt_flag
	alpha2 = args.alpha2
	beta2 = args.beta2
	pval = args.pvalue
	chrM_conta = args.chrM_contaminant

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
	META_DICT, ALL_CELL_TYPES = meta_to_dict(meta_file,tissue)
	DICT_VARIANTS = build_dict_variants(variant_file,window)
	
	# 3. Code to run in parallel all bins
	if (CORE > 1):
		pool = mp.Pool(CORE)
		
		# Step 3.1: Use loop to parallelize
		for row in DICT_VARIANTS.keys():
			# This funtion writes in temp files the results
			pool.apply_async(run_interval, args=(row,DICT_VARIANTS,META_DICT,BAM,FASTA,tmp_dir,BQ,MQ,ALT_FLAG,alpha2,beta2,pval,chrM_conta), callback=collect_result)
				   
		# Step 3.2: Close Pool and let all the processes complete	
		pool.close()
		pool.join()
	else:
		for row in DICT_VARIANTS.keys():
			# This funtion writes in temp files the results
			collect_result(run_interval(row,DICT_VARIANTS,META_DICT,BAM,FASTA,tmp_dir,BQ,MQ,ALT_FLAG,alpha2,beta2,pval,chrM_conta))
	
	# 4. Write "Long" file
	concatenate_sort_temp_files_and_write(out_prefix, tmp_dir)

	# 4.5 Compute fusion matrix
	if fusion_file:
		fusions = collect_cells_with_fusions(fusion_file)
	else:
		fusions = pd.DataFrame()

	# 5. Write matrixes ("wide") DP/ALT/VAF/Binarized files
	pivot_long_dataframe(out_prefix,fusions)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



