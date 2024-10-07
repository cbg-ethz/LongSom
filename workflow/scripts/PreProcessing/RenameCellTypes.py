import pandas as pd
import argparse

def rename_cell_types(input,output,cancer_cell_type):
    barcodes = pd.read_table(input)
    barcodes['Input_cell_type'] = barcodes['Cell_type']
    barcodes['Cell_type'] = barcodes['Cell_type'].apply(lambda x: 'Cancer' if x==cancer_cell_type else 'Non-Cancer')
    barcodes.to_csv(output,sep='\t',index=False)


def initialize_parser():
    parser = argparse.ArgumentParser(description='Rename cell types in cancer/non-cancer')
    parser.add_argument('--input', type=str, default=1, help='User input barcode file', required = True)
    parser.add_argument('--output', type=str, default=1, help='Barcode file with redefined cancer/non-cancer celltypes', required = True)
    parser.add_argument('--cancer_cell_type', type=str, default = 'Sample', help='User defined cancer cell type', required = False)
    return (parser)


def main():
    # 1. Arguments
    parser = initialize_parser()
    args = parser.parse_args()

    input=args.input
    output=args.output
    cancer_cell_type=args.cancer_cell_type

    # 2. Rename cell types
    rename_cell_types(input,output,cancer_cell_type)


#-------------
# Execute code
#-------------

if __name__ == '__main__':
    main()


