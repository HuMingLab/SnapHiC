
import pandas as pd
import numpy as np
import scipy as sp
import networkx as nx
import os
import re
import sys
import random
from collections import deque
import gc


CHROM_DICT = {"chr19" : 61431566}
SETNAME=["MY_214"]
#START = 30e6
#END = 32e6
BIN = 1e4
ALPHA = 0.05
DIST = 1e6
edge_filename = "/Users/abnousa/data/single_cell/snap_hic/data/inputs/miao_mESC_unphased/binned_data"


# In[3]:

#@profile
def get_rwr(edge_filename, binsize = BIN, distance = DIST, chrom = list(CHROM_DICT.keys())[0], chrom_len = CHROM_DICT[list(CHROM_DICT.keys())[0]], alpha = ALPHA, final_try = False):
    gc.collect()
    edgelist = pd.read_csv(edge_filename, sep = "\t", header = None, names = ["chr1", "x1", "x2", "chr2", "y1", "y2"])
    edgelist = edgelist[(edgelist['chr1'] == chrom) & (edgelist['chr1'] == edgelist['chr2'])]
    edgelist = bin_matrix(edgelist, binsize)
    NUM = int(np.ceil(chrom_len / binsize))
    
    edges = pd.DataFrame({'x1':list(range(1, NUM)), 'y1':list(range(2, NUM+1))})
    edges = pd.concat([edges, edgelist[['x1', 'y1']]], axis = 0)
    edges.loc[:,'weight'] = 1
    
    g = get_stochastic_matrix_from_edgelist(edges)
    r = solve_rwr(g, alpha, final_try)
    if r.isinstance(str) and r == "try_later":
        return r

    df = reformat_sparse_matrix(r, binsize, distance)
    df.loc[:,'chr1'] = chrom
    df.loc[:,'chr2'] = chrom
    df = df[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'value']]
    return df


def bin_matrix(df, binsize):
    df.loc[:,'x1'] = df.loc[:,'x1'] // binsize
    df.loc[:,'y1'] = df.loc[:,'y1'] // binsize
    return df


def get_stochastic_matrix_from_edgelist(edgelist):
    gc.collect()
    g = nx.from_pandas_edgelist(edgelist, source = 'x1', target = 'y1', edge_attr = ['weight'], create_using = nx.Graph())
    degrees = np.array([g.degree(i) for i in g.nodes()])
    m = sp.sparse.csc_matrix(nx.adjacency_matrix(g).astype(float))
    m.data /= degrees[m.indices] #stochastic matrix
    return m


def solve_rwr(stoch_matrix, alpha = ALPHA, final_try = False):
    gc.collect()
    m = stoch_matrix*(1-alpha)
    m = m.transpose()
    y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")
    A = y - m
    try:
        s = sp.sparse.linalg.spsolve(A, y)
    except:
        A = A.todense()
        y = y.todense()
        s = sp.linalg.solve(A, y)
        s = sp.sparse.csr_matrix(s)
    finally:
        if final_try:
            raise Exception("Cannot allocate enough memory of solving RWR")
        else:
            return "try_later"
    s *= alpha
    s += s.transpose()
    return s.tocoo()


def reformat_sparse_matrix(m, binsize, distance):
    max_bin_distance = distance // binsize
    df = pd.DataFrame({'x1':m.row, 'y1':m.col, 'value':m.data})
    df = df[((df['y1'] - df['x1']) <= max_bin_distance) & ((df['y1'] - df['x1']) > 0)]
    df.iloc[:,0] = df.iloc[:,0] * binsize
    df.iloc[:,1] = df.iloc[:,1] * binsize
    df.loc[:,'x2'] = df['x1'] + binsize
    df.loc[:,'y2'] = df['y1'] + binsize
    df[['x1', 'x2', 'y1', 'y2']] = df[['x1', 'x2', 'y1', 'y2']].astype(int) 
    return df

def normalize_along_diagonal(d, trim = 0.01):
    trim_value = d['value'].quantile(1-trim)
    #print(trim_value)
    remaining = d[d['value'] < trim_value]['value']
    mu = np.mean(remaining)
    sd = np.std(remaining)
    #print(trim_value, remaining.shape[0], mu, sd)
    d.loc[:,'value'] = (d['value'] - mu) / sd
    return d

def determine_proc_share(indir, chrom_lens, n_proc, rank):
    filenames = os.listdir(indir)
    filenames = [name for name in filenames if name.endswith(".bedpe")]
    filenames.sort()
    #print(filenames)
    setnames = [re.search(r".*?\.", fname).group()[:-1] for fname in filenames]
    chrom_names = list(chrom_lens.keys())
    chrom_names.sort()
    jobs = [(chrom_names[i], filenames[j], setnames[j]) for i in range(len(chrom_names)) for j in range(len(filenames))]
    
    indices = list(range(rank, len(jobs), n_proc))
    proc_jobs = [jobs[i] for i in indices]
    #proc_jobs = deque(proc_jobs)
    #proc_jobs.rotate(rank)
    #proc_jobs = list(proc_jobs)
    random.shuffle(proc_jobs)
    return proc_jobs
    

def get_rwr_for_all(indir, outdir = None, binsize = BIN, alpha = ALPHA, dist = DIST, chrom_lens = CHROM_DICT, normalize = False, n_proc = 1, rank = 0, genome = 'mouse'):
    if not outdir:
        outdir = indir
        outdir = os.path.join(outdir, "rwr")
    try:
        os.makedirs(outdir)
    except:
        pass

    processor_jobs = determine_proc_share(indir, chrom_lens, n_proc, rank)
    retry_filename = os.path.join(outdir, ('_'.join([str(rank), "retry", "instances"]) + ".txt"))
    attempt_counter = 0
    attempts_allowed = 10
    while len(processor_jobs > 0):
        for chrom, filename, setname in processor_jobs:
            filepath = os.path.join(indir, filename)
            final_try = False if attempt_counter < attempts_allowed else True
            df = get_rwr(filepath, binsize = binsize, distance = dist, chrom = chrom, chrom_len = chrom_lens[chrom], alpha = alpha)
            if isinstance(df, str) and df == "try_later":
                with open(retry_filename, 'a') as ofile:
                    ofile.write("\t".join([chrom, filename, setname]) + "\n")
                continue
            output_filename = os.path.join(outdir, ".".join([setname, chrom, "rwr", "bedpe"]))
            df.sort_values(['x1', 'y1'], inplace = True)
            df.to_csv(output_filename, sep = "\t", header = None, index = False)
            if normalize:
                output_filename = os.path.join(outdir, ".".join([setname, chrom, "normalized", "rwr", "bedpe"]))
                df = df.groupby(df['y1'] - df['x1'], as_index = False).apply(normalize_along_diagonal).reset_index(drop = True)
                df.sort_values(['x1', 'y1'], inplace = True)
                df.to_csv(output_filename, sep = "\t", header = None, index = False)
        try:
            print("rank", rank, ": attempting to rerun failed jobs. Attempt #", attempt_counter + 1)
            sys.stdout.flush()
            with open(retry_filename, 'r') as infile:
                jobs = infile.readlines()
            jobs = [line.split() for line in jobs]
            processor_jobs = jobs
            os.remove(retry_filename)
            attempt_counter += 1
        except:
            print("rank", rank, ": no remaining jobs or parser failed")
            sys.stdout.flush()

if __name__ == "__main__":
    get_rwr_for_all(INDIR, normalize = True)

