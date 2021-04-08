import h5py
import numpy as np
import pandas as pd
import os
import sys
import tqdm
from tqdm import trange
import argparse
from scipy.sparse import coo_matrix

def main(args):
    ## hdf5 files should be in higashi format: one dataset called "coordinates" for coordinates plus one dataset celled "cell_x" for each cell "x" (x is an integer)
    ## bedpe files are upper triangular. i.e. the first coordiante is smaller than the second coordinate on each line and 
    ## bedpe files should include a value for the very last binpair (even if it's 0). The matrix size is determined based on the largest index in the bedpe file.
    if args.filetype not in ['npz', 'hdf5', 'bedpe']:
        raise Exception("Filetype should be npz, hdf5(higashi), or bedpe")
    if filetype == "hdf5":
        process_hdf_file(args.input_filename, args.outdir, args.chrom, binsize = args.binsize, dist = args.dist)

def read_imputed_matrix(filename):
    with h5py.File(filename, "r") as impute_f:
        coordinates = impute_f['coordinates']
        xs, ys = coordinates[:, 0], coordinates[:, 1]
        size = int(np.max(ys)) + 1
        cell_list = trange(len(list(impute_f.keys())) - 1)
        print("cells in file:", len(cell_list))
        m1 = np.zeros((size, size))
        for i in cell_list:
            m1 *= 0.0
            proba = np.array(impute_f["cell_%d" % i])
            m1[xs.astype('int'), ys.astype('int')] += proba
            m1 = m1 + m1.T	
            yield i, m1

def normalize_along_diagonal_from_numpy(d, chrom, max_bin_distance, output_filename, binsize, remove_bins, trim = 0.01):
    #df_all = pd.DataFrame()
    if os.path.exists(output_filename):
        os.remove(output_filename)
    with open(output_filename, "a") as f:
        for offset in range(1, max_bin_distance + 1):
            r, c = get_nth_diag_indices(d, offset)
            vals_orig = d[r,c].tolist()

            if isinstance(vals_orig[0], list):
                vals_orig = vals_orig[0]
            #print(type(vals_orig))
            #print(type(vals_orig[0]))
            #print(len(vals_orig))
            #vals_orig = vals_orig[0]

            vals = vals_orig.copy()
            vals.sort()
            vals.reverse()
            trim_value = vals[(round(trim * len(vals)) - 1)]
            trim_index = round(trim * len(vals)) - 1
            remaining = vals[(trim_index):]
            mu = np.mean(remaining)
            sd = np.std(remaining)
            #print(vals_orig[:5])
            #print('musd')
            #print(mu, sd)
            #vals_orig = vals_orig[0][:]
            vals_orig = (np.array(vals_orig) - mu) / sd
            if sd < 1e-6:
                vals_orig = [0] * len(vals_orig)
                #print(len(vals_orig))
            #print(offset,mu,sd,max(vals_orig))
            #print(vals_orig[:5])
            #print(r.shape, c.shape, vals_orig.shape)
            df = pd.DataFrame({'x1': r, 'y1': c, 'v': vals_orig})
            df['x1'] = ((df['x1'] + remove_bins) * binsize).astype(int)
            df['y1'] = ((df['y1'] + remove_bins) * binsize).astype(int)
            df['x2'] = (df['x1'] + binsize).astype(int)
            df['y2'] = (df['y1'] + binsize).astype(int)
            df['chr1'] = chrom
            df['chr2'] = chrom
            df = df[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'v']]
            df.to_csv(f, mode='a', header=False, sep = "\t", index = False)

def get_nth_diag_indices(mat, offset):
    rows, cols_orig = np.diag_indices_from(mat)
    cols = cols_orig.copy()
    if offset > 0:
        cols += offset
        rows = rows[:-offset]
        cols = cols[:-offset]
    return rows, cols

def process_hdf_file(input_filename, outdir, chrom, binsize = 1e4, dist = 2e6):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cells = read_imputed_matrix(input_filename)
    for cell_name, mat in cells:
        output_filename = f"{outdir}/{cell_name}.{chrom}.normalized.rwr.bedpe"
        buffer_size = 50000
        max_bin_distance = int((dist + buffer_size) // binsize)
        remove_bins = int(buffer_size // binsize)
        normalize_along_diagonal_from_numpy(mat, chrom, max_bin_distance, output_filename, binsize, remove_bins)

def process_npz_files(indir, outdir, chrom, binsize = 1e4, dist = 2e6):
    if not os.path.exists(outdir):
        os.makefirs(outdir)
    cells = glob.glob(f"{indir}/*{chrom}*{.npz}")
    for cell_number, filename in enumerate(cells):
        mat = np.load(filename)
        output_filename = f"{outdir}/cell{cell_number}.{chrom}.nomralized.rwr.bedpe"
        buffer_size = 50000
        max_bin_distance = int((dist + buffer_size) // binsize)
        remove_bins = int(buffer_size // binsize)
        normalize_along_diagonal_from_numpy(mat, chrom, max_bin_distance, output_filename, binsize, remove_bins)

def process_bedpe_files(indir, outdir, chrom, binsize = 1e4, dist = 2e6):
    if not os.path.exists(outdir):
        os.makefirs(outdir)
    cells = glob.glob(f"{indir}/*{chrom}*{.bedpe}")
    for cell_number, filename in enumerate(cells):
        df = pd.read_csv(filename, sep = "\t", header = None)
        mat = convert_bedpe_to_matrix(df, binsize)
        output_filename = f"{outdir}/cell{cell_number}.{chrom}.nomrlalized.rwr.bedpe"
        buffer_size = 50000
        max_bin_distance = int((dist + buffer_size) // binsize)
        remove_bins = int(buffer_size // binsize)
        normalize_along_diagonal_from_numpy(mat, chrom, max_bin_distance, output_filename, binsize, remove_bins)

def convert_bedpe_to_matrix(df, binsize):
    df.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'v']
    df['x1'] = (df['x1'] // binsize).astype(int)
    df['y1'] = (df['y1'] // binsize).astype(int)
    m = coo_matrix((df['v'], (df['x1'], df['y1'])), shape = (df['y1'].max() + 1, df['y1'].max() + 1))
    m = m.todense()
    return m

def setup_parser():
    parser = argparse.ArgumentParser(description='Normalizing imputed snHiC data')
    parser.add_argument('-i', '--infile', help='file computing imputed data')
    parser.add_argument('-t', '--filetype', help = "input filetype; should be one of: npz, hdf5, bedpe")
    parser.add_argument('-o', '--outdir', help = "output directory for output files")
    parser.add_argument('-b', '--binsize', help = "binsize used for imputation. default 1e04."), default = 1e4, type = int)
    parser.add_argument('-d', '--dist', help = "maximum distance between a pair of bins", default = 2e6)
    parser.add_argument('-c', '--chrom', help = "chromosome to process", required = True)
    return parser

args = parser.parse_args()

if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args()
    main(args)
