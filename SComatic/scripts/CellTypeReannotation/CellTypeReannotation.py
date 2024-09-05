import timeit
import argparse
import pandas as pd
from collections import Counter

def collect_cells_with_SNVs(SNV_file, min_variants):
	SNVs = pd.read_csv(SNV_file, sep='\t')
	SNVs['INDEX'] = SNVs['#CHROM'] + ':' + SNVs['Start'].astype(str) + ':' + SNVs['ALT_expected'].str.split(',', n=1, expand=True)[0]
	
	# Keeping only cells with at least min_variants sites covered
	BC_COV = SNVs[SNVs['VAF'] != '.']
	BC_COV = Counter(BC_COV['CB'])
	BC_COV_min = [k for k,v in BC_COV.items() if v>=min_variants]
	SNVs = SNVs[SNVs['CB'].isin(BC_COV_min)]

	# Filter only cells with called mutations
	SNVs = SNVs[SNVs['MutationStatus'] == 'PASS']
	
	# Collect mutated barcodes
	return list(SNVs['CB']),BC_COV,BC_COV_min

def collect_cells_with_fusions(fusion_file):
	fusions = pd.read_csv(fusion_file, sep='\t')

	# Create Fusion:barcode index to remove duplicates
	fusions['INDEX'] = fusions['FusionName'] + ':' + fusions['barcodes']

	# Drop duplicates
	fusions = fusions.drop_duplicates(subset='INDEX', keep="last")

	# Collect mutated barcodes
	return list(fusions['barcodes'])


def collect_cancer_cells(cells_with_SNVs,cells_with_fusions,BC_COV,min_variants,min_frac):

	cells = cells_with_SNVs + cells_with_fusions

	# Count variants per cell
	variants_per_cell = Counter(cells)
	# Compute the fraction of variants with coverage that is mutated in each cell
	mean_variants_mut_per_cell = { k : (v/BC_COV[k] if BC_COV[k]>=min_variants else 0) for k,v in variants_per_cell.items()}

	# Determine cancer cells base on user-defined # of variants and fraction of mutated variants
	#cancer_cells = [k for k,v in variants_per_cell.items() if v>=min_variants if mean_variants_mut_per_cell[k]>=min_frac]
	cancer_cells = [k for k,v in variants_per_cell.items() if mean_variants_mut_per_cell[k]>=min_frac]

	return cancer_cells


def write_reannotated_cell_types(cancer_cells,BC_COV_min,bc_file,out_file):
	bcs = pd.read_csv(bc_file, sep='\t')
	bcs = bcs[bcs['Index'].isin(BC_COV_min)]
	
    # Save automated annotation
	bcs['Automated_Cell_type'] = bcs['Cell_type']

    # Cell reannotation based on HCCVs
	bcs['Reannotated_cell_type'] = ['Cancer' if i in cancer_cells else 'NonCancer' for i in bcs['Index']]
	
    # Copy reannotated celltypes for SplitBamCellTypes.py compatibility
	bcs['Cell_type'] = bcs['Reannotated_cell_type']

	bcs['Cancer_Color'] = bcs['Reannotated_cell_type'].apply(lambda x: '#94C773' if x=='NonCancer' else '#8F79A1')
	
    # Writing output file
	bcs.to_csv(out_file, sep='\t', index = False)


def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to get the alleles observed in each unique cell for the variant sites')
	parser.add_argument('--SNVs', type=str, help='HCCV SNV calls (obtained by SingleCellGenotype.py), tsv file', required = True)
	parser.add_argument('--fusions', type=str, help='HCCV fusion calls (obtained by CTAT_Fusion.smk), tsv file', required = True)
	parser.add_argument('--outfile', type=str, help='Output tsv file', required = True)
	parser.add_argument('--meta', type=str, help='Barcodes tsv file', required = True)
	parser.add_argument('--min_variants', type=int, default = 3, help='Minimum # of variants covered to consider a cell', required = False)
	parser.add_argument('--min_frac', type=float, default = 0.2, help='Minimum fraction of covered variants being mutated to call a cancer cell', required = False)

	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	SNV_file = args.SNVs
	fusion_file = args.fusions
	out_file = args.outfile
	bc_file = args.meta
	min_variants = args.min_variants
	min_frac = args.min_frac

	# Set outfile name
	print("Outfile: " , out_file ,  "\n") 

	# 1. Collect cell barcodes where at least a SNV HCCV was found 
	cells_with_SNVs,BC_COV,BC_COV_min = collect_cells_with_SNVs(SNV_file,min_variants)

	# 2. Collect cell barcodes where at least a fusion HCCV was found
	cells_with_fusions = collect_cells_with_fusions(fusion_file)
	
	# 3. Find cancer cells (cells with the user defined min. # of HCCV variants)
	cancer_cells = collect_cancer_cells(cells_with_SNVs,cells_with_fusions,BC_COV,min_variants,min_frac)

	# 4. Write reannotated celltype file
	write_reannotated_cell_types(cancer_cells,BC_COV_min,bc_file,out_file)


if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



