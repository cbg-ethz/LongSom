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
	fi=pysam.Samfile(args.input, "rb",threads=args.cpu)

	with open(args.barcodes,'r') as bc_file:
		barcodes = [i[:-1] for i in bc_file.readlines()]

	bcout_dict = {}

	for bc in barcodes:
		bcout = pysam.Samfile('barcodes/'+args.sample+"/"+bc+".bam", "wb", template = fi)
		bcout_dict[bc] = bcout

	for READ in fi.fetch():
		BC = READ.get_tag('CB')
		try:
			bcout_dict[BC].write(READ)
		except KeyError:
			continue

	for bc,bcout in bcout_dict.items():
		bcout.close()
		pysam.index('barcodes/'+args.sample+"/"+bc+".bam")

	with open('barcodes/'+args.sample+'/terminado_splitbc.txt','w') as f:
		f.write('Tout le monde descend')

def parse_args():
	parser = argparse.ArgumentParser(
		prog='filter_dist_nontum_final.py', 
		usage='python3 FMI_mutations_to_cells.py --bam <sorted.bam> --csv FMI_mutations.csv',
		description='Creates csv files of mutated cell names and their mutations'
	)

	parser.add_argument(
			'--cpu', type=int, default=8,
			help='number of cpus for pysam'
	)
	parser.add_argument(
			'--sample', type=str,
			help='sample name'
	)
	parser.add_argument(
			'--barcodes', type=str,
			help='barcodes file'
	)
	parser.add_argument(
			'--input', type=str,
			help='bam input'
	)

	args = parser.parse_args()
	return args

if __name__ == '__main__':
	args = parse_args()
	main(args)
