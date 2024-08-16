import timeit
import argparse
import pandas as pd
from collections import Counter

def collect_cells_with_SNVs(SNV_file):
	SNVs = pd.read_csv(SNV_file, sep='\t')
	BC_COV= SNVs[SNVs['VAF'] != '.']['CB'].unique()
	SNVs = SNVs[SNVs['CB'].isin(BC_COV)]
	# Filter only cells with called mutations
	SNVs = SNVs[SNVs['MutationStatus'] == 'PASS']
	

	# Collect mutated barcodes
	return list(SNVs['CB']),BC_COV

def collect_cells_with_fusions(fusion_file):
	fusions = pd.read_csv(fusion_file, sep='\t')

	# Create Fusion:barcode index to remove duplicates
	fusions['INDEX'] = fusions['FusionName'] + ':' + fusions['barcodes']

	# Drop duplicates
	fusions = fusions.drop_duplicates(subset='INDEX', keep="last")

	# Collect mutated barcodes
	return list(fusions['barcodes'])


def collect_cancer_cells(cells_with_SNVs,cells_with_fusions,min_variants):
	cells = cells_with_SNVs + cells_with_fusions

	# Count variants per cell
	variants_per_cell = Counter(cells)

	# Determine cancer cells base on user-defined # of variants
	cancer_cells = [k for k,v in variants_per_cell.items() if v>=min_variants]
	rescue_cells = [k for k,v in variants_per_cell.items() if v>=(min_variants-1)]

	return cancer_cells,rescue_cells


def write_reannotated_cell_types(cancer_cells,rescue_cells,BC_COV,bc_file,out_file):
	bcs = pd.read_csv(bc_file, sep='\t')
	bcs = bcs[bcs['Index'].isin(BC_COV)]
	meta_dict = dict(zip(bcs['Index'],bcs['Cell_type']))
	rescue_cells = [i for i in bcs['Index'] if i in rescue_cells]
	rescue_cells = [i for i in rescue_cells if meta_dict[i]=='HGSOC']
	cancer_cells = cancer_cells+rescue_cells
	
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
	parser.add_argument('--min_variants', type=int, default = 2, help='Minimum # of variants to call a cancer cell', required = False)

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

	# Set outfile name
	print("Outfile: " , out_file ,  "\n") 

	# 1. Collect cell barcodes where at least a SNV HCCV was found 
	cells_with_SNVs,BC_COV = collect_cells_with_SNVs(SNV_file)

	# 2. Collect cell barcodes where at least a fusion HCCV was found
	cells_with_fusions = collect_cells_with_fusions(fusion_file)
	
	# 3. Find cancer cells (cells with the user defined min. # of HCCV variants)
	cancer_cells,rescue_cells = collect_cancer_cells(cells_with_SNVs,cells_with_fusions,min_variants)

	# 4. Write reannotated celltype file
	write_reannotated_cell_types(cancer_cells,rescue_cells,BC_COV,bc_file,out_file)


if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')



