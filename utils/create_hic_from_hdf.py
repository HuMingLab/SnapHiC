import h5py
import numpy as np
import pandas as pd
import glob
import os

chroms=['chr' + str(i) for i in range(1,20)]
rwr_dir = "/oasis/tscc/scratch/abnousa/snapHiC/outputs/fraser_binsizes/fraser_5kb_subsamples/fraser_700/rwr"
hic_dir = "/oasis/tscc/scratch/abnousa/snapHiC/outputs/fraser_binsizes/fraser_5kb_subsamples/fraser_700/hic"
num_cells = 700 

hic_input_fname = f"{hic_dir}/allChr.hic.input"
if os.path.exists(hic_input_fname):
	os.remove(hic_input_fname)
print(hic_input_fname)
for chrom in chroms:
	print(chrom)
	hf = h5py.File(f'{hic_dir}/{chrom}.normalized.combined.bedpe.cells.hdf', 'r')
	d = hf[chrom]
	d = d[:,:]
	outlier_counts = np.apply_along_axis(arr = d, axis = 1, func1d = lambda x: len(np.where(x>1.96)[0]))
	d = pd.read_csv(glob.glob(f'{rwr_dir}/*{chrom}.normalized.rwr.bedpe')[0], sep = "\t", header = None)
	d[6] = np.ceil(outlier_counts/num_cells * 100).astype(int)
	for i in [1,2,4,5,6]:
		d[i] = d[i].astype(int)
	d.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'score']
	d['str1'] = 0
	d['str2'] = 1
	d['frag1'] = 0
	d['frag2'] = 1
	d = d[['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score']]
	d.to_csv(hic_input_fname, sep = "\t", header = False, index = False, mode = 'a')
