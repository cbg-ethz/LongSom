import pandas as pd
import argparse

def main(args):
    filtered = pd.read_csv(args.filtered, names=['barcodes'], sep='\t')
    raw = pd.read_csv(args.raw, names=['barcodes'],compression='gzip', encoding = "ISO-8859-1", sep='\t')
    filtered_set = set(filtered['barcodes'])
    raw_list = list(raw['barcodes'])


    emptydrops = {cell[:-2] for cell in raw_list if cell not in filtered_set}

    csv = pd.DataFrame({'barcodes': list(emptydrops)})
    csv.to_csv('emptydroplets/barcodes/{}.tsv'.format(args.sample), sep='\t', index=False)

def parse_args():
    parser = argparse.ArgumentParser(
        prog='split_cells_bam.py', 
        usage='python3 split_cells_bam.py --bc_dir <bc_dir> --sample <sample> ',
        description='Divides bamfiles per cell prior to UMI deduplication'
    )
    parser.add_argument(
        '--filtered', type=str,
        help='Absolute or relative path(s) to directory containing barcodes files sample.whatever.txt, tsv format'
    )
    parser.add_argument(
        '--raw', type=str,
        help='sample name (should not contain ".")'
    )
    parser.add_argument(
        '--sample', type=str,
        help='sample name (should not contain ".")'
    )
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)



