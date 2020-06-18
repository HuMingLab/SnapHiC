
import os
import numpy as np
import pandas as pd
from itertools import repeat
import re


NUM_PROCESSOR = 1
RANK = 0
BIN = 1e4
LOW_CUTOFF = 5e3
INDIR = "/Users/abnousa/data/single_cell/snap_hic/data/inputs/miao_mESC_unphased"
FILE_SUFFIX = "valid_pairs.rm_hotspot.sorted.txt"
CHR_COLUMNS = [2, 6]
POS_COLUMNS = [3, 7]


def extract_setname(filepath, suffix):
    return re.sub("\.*" + suffix, "", os.path.basename(filepath))


def get_filepaths(input_dir, suffix = "valid_pairs.rm_hotspot.sorted.txt"):
    filenames = [os.path.join(input_dir, name) for name in os.listdir(input_dir) if name.endswith(suffix)]
    return filenames


def get_proc_filenames(filenames, n_proc = NUM_PROCESSOR, rank = RANK):
    indices = list(range(rank, len(filenames), n_proc))
    proc_fnames = [filenames[i] for i in indices]
    return proc_fnames


def bin_file(filename, binsize, outdir, chr_columns, pos_columns, file_suffix, low_cutoff):
    df = pd.read_csv(filename, sep = "\t", header = None)   
    
    chr_columns = [i-1 for i in chr_columns]
    pos_columns = [i-1 for i in pos_columns]
    df = df.iloc[:, chr_columns + pos_columns]
    
    df.columns = ["chr1", "chr2", "x1", "y1"]
    #print(df.head())
    
    #remove low cutoff distance read pairs
    df = df[abs(df['x1'] - df['y1']) >= low_cutoff]
    df['x1'] = (df['x1'] // binsize * binsize).astype(int)
    df['y1'] = (df['y1'] // binsize * binsize).astype(int)
        
    #keep intra chr
    df = df[(df['chr1'] == df['chr2'])]

    #keep autosomal
    df = df[df['chr1'].isin(["chr" + str(i) for i in range(23)])]
    
    #make x1 be the smaller bin
    df.loc[df['x1'] > df['y1'], ['x1', 'y1']] = df.loc[df['x1'] > df['y1'], ['y1', 'x1']].values
    
    #remove duplicate rows
    df.drop_duplicates(inplace = True)
    
    #remove adjacent binpairs
    df = df[df['y1'] - df['x1'] >= 2 * binsize]
    
    df['x2'] = (df['x1'] + binsize).astype(int)
    df['y2'] = (df['y1'] + binsize).astype(int)
    df = df[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']]
    
    name = extract_setname(filename, file_suffix)
    
    df.to_csv(os.path.join(outdir, name + ".bedpe"), sep = "\t", header = None, index = False)
    return name

def bin_sets(indir, file_suffix, binsize = 1e4, outdir = None, chr_columns = [2, 6], pos_columns = [3, 7], low_cutoff = 5e3, n_proc = 1, rank = 0):
    if not outdir:
        outdir = indir
        outdir = os.path.join(outdir, "binned_data")
    filenames = get_filepaths(indir, file_suffix)
    proc_filenames = get_proc_filenames(filenames, n_proc, rank)
    try:
        os.makedirs(outdir)
    except:
        pass

    #setnames = []
    for filename in proc_filenames:
        setname = bin_file(filename, binsize, outdir, chr_columns, pos_columns, file_suffix, low_cutoff)
    #    setnames.append(setname)
    #return setnames

if __name__ == "__main__":
    bin_sets(INDIR, FILE_SUFFIX)
