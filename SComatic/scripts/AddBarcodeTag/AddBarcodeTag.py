import pysam
import argparse

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq)
    letters = [complement[base] for base in bases] 
    letters = ''.join(letters)
    reverse = letters[::-1]
    return reverse

def main(args):
    fi=pysam.AlignmentFile(args.input, "rb", threads=args.cpu)
    fo=pysam.AlignmentFile(args.output, "wb", threads=args.cpu, template=fi)

    for READ in fi.fetch():
        NAME = READ.query_name
        BC = NAME.split('.')[-1]
        BC = reverse_complement(BC)
        READ.tags += [('CB', BC)]
        fo.write(READ)

    fo.close()
    fi.close()


def parse_args():
    parser = argparse.ArgumentParser(
        prog='AddBarcodeTag.py', 
        usage='python3 ../../bin/SComatic/scripts/AddBarcodeTag.py --input <in.bam> --output <out.bam>',
        description='Adds barcode tag to bam file, from read name'
    )

    parser.add_argument(
            '--cpu', type=int, default=8,
            help='number of cpus for pysam'
    )
    parser.add_argument(
            '--input', type=str,
            help='bam input'
    )
    parser.add_argument(
            '--output', type=str,
            help='bam output'
    )

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
