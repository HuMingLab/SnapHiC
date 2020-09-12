import pandas as pd
import numpy as np
import argparse
import sys

def main():
	parser = create_parser()
	args = parser.parse_args()
	output_filename=f'{args.indir}/hic/{args.setname}.normalized.hic'
	filename = f'{args.indir}/interactions/combined_significances.bedpe'
	d = pd.read_csv(filename, sep = "\t")
	d['neglog10_fdr'] = -np.log10(d['fdr_dist'] + 1e-6)
	#print(d.describe())
	d['score'] = round(d['neglog10_fdr']).astype(int)
	d['str1'] = 0
	d['str2'] = 1
	d['frag1'] = 0
	d['frag2'] = 1
	#print(d['neglog10_fdr'].describe())
	d = d[['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score']]
	#print(d['score'].describe())
	sys.stdout.flush()
	d.to_csv(output_filename + ".input", index = False, header = None, sep = "\t")

def create_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--indir', action = 'store', required = True, \
        	                help = 'input directory; should be same is outdir set for snapHiC')
	parser.add_argument('-s', '--setname', action = 'store', required = True, \
                                help = 'name that will be used for naming the output file')
	return parser

main()
