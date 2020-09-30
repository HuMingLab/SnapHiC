
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
#import dask
#import dask.array
#import scikits.umfpack


CHROM_DICT = {"chr19" : 61431566}
SETNAME=["MY_214"]
#START = 30e6
#END = 32e6
BIN = 1e4
ALPHA = 0.05
DIST = 1e6
edge_filename = "/Users/abnousa/data/single_cell/snap_hic/data/outputs/miao_mESC_unphased/binned"


# In[3]:

#@profile
def get_rwr_en(edge_filename, binsize = BIN, distance = DIST, chrom = list(CHROM_DICT.keys())[0], chrom_len = CHROM_DICT[list(CHROM_DICT.keys())[0]], \
                                                alpha = ALPHA, final_try = False, logfile = None, parallel = False, threaded_lock = None):
    gc.collect()
    setname = edge_filename[(edge_filename.rfind('/') + 1):]
    #print('computing rwr for', setname, chrom)
    #sys.stdout.flush()
    edgelist = pd.read_csv(edge_filename, sep = "\t", header = None, names = ["chr1", "x1", "x2", "chr2", "y1", "y2"])
    edgelist = edgelist[(edgelist['chr1'] == chrom) & (edgelist['chr1'] == edgelist['chr2'])]
    edgelist = bin_matrix(edgelist, binsize)
    NUM = int(np.ceil(chrom_len / binsize))
    #print('NUM', NUM)
    
    edges = pd.DataFrame({'x1':list(range(0, NUM-1)), 'y1':list(range(1, NUM))})
    edges = pd.concat([edges, edgelist[['x1', 'y1']]], axis = 0)
    edges.loc[:,'weight'] = 1
    #print('deleting edgelist', setname, chrom)
    #sys.stdout.flush()
    del edgelist
    #print('getting stoc matrix', setname, chrom)
    #sys.stdout.flush()
    g = get_stochastic_matrix_from_edgelist(edges)
    #print('matrix returned', setname, chrom)
    #sys.stdout.flush()
    gc.collect()
    #print('solving rwr', setname, chrom)
    #sys.stdout.flush()
    msg = f'solved equation for {chrom} {setname}\n'
    r = solve_rwr_en(g, alpha, final_try, setname, chrom)
    if logfile and parallel:
        logfile.Write_shared(msg.encode('utf-8'))
        logfile.Sync()
    elif logfile:
        if threaded_lock:
            threaded_lock.acquire()
        logfile.write(msg)
        logfile.flush()
        if threaded_lock:
            threaded_lock.release()
    #print('r returned', setname, chrom)
    #print(type(r), setname, chrom)
    #sys.stdout.flush()
    del g, edges
    gc.collect()
    if isinstance(r, str) and r == "try_later":
        #print('returning try later', setname, chrom)
        return r
    #print('reformatting', setname, chrom)
    return r


def bin_matrix(df, binsize):
    df.loc[:,'x1'] = df.loc[:,'x1'] // binsize
    df.loc[:,'y1'] = df.loc[:,'y1'] // binsize
    return df


def get_stochastic_matrix_from_edgelist(edgelist):
    g = nx.from_pandas_edgelist(edgelist, source = 'x1', target = 'y1', edge_attr = ['weight'], create_using = nx.Graph())
    degrees = np.array([g.degree(i) for i in g.nodes()])
    m = sp.sparse.csc_matrix(nx.adjacency_matrix(g).astype(float))
    m.data /= degrees[m.indices] #stochastic matrix
    del g, degrees
    return m


def solve_rwr_en(stoch_matrix, alpha = ALPHA, final_try = False, setname = None, chrom = None):
    gc.collect()
    #print('first', setname, chrom); sys.stdout.flush()
    m = stoch_matrix*(1-alpha)
    #print('second', setname, chrom)
    #sys.stdout.flush()
    m = m.transpose()
    #print('thhird', setname, chrom); sys.stdout.flush()
    y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")
    #print('fifth', setname, chrom); sys.stdout.flush()
    A = y - m
    #print('trying', setname, chrom)
    #sys.stdout.flush()
    #try:
    #    s = None
    #    #s = sp.sparse.linalg.spsolve(A, y)
    #    s = scikits.umfpack.spsolve(A, y)
    #    if isinstance(s, np.ndarray):
    #        s = sp.sparse.csr_matrix(s)
    #    #print('s first try', setname, chrom)
    #    #sys.stdout.flush()
    #except:
    #    #print('attempting to delete', setname, chrom)
    #    #sys.stdout.flush()
    #    if s is not None:
    #        del s
    #    #print("sparse solver failed. trying dense solver.")
    #    #sys.stdout.flush()
    #    gc.collect()
    try:
        s = None
        #print('second try', setname, chrom)
        #sys.stdout.flush()
        #A = A.todense()
        #A = dask.array.from_array(A.todense())
        A = A.todense()
        #print('have a', setname, chrom)
        y = y.todense()
        #y = dask.array.from_array(y.todense())
        #print('have y', setname, chrom, A.shape, y.shape)
        s = sp.linalg.solve(A, y)
        #s = dask.array.linalg.solve(A, y)
        #print('working on s', setname, chrom)
        #print('type1', type(s))
        #s = s.compute()
        #print('type2', type(s))
        #s = sp.sparse.csr_matrix(s)
        #print('s on second try;', setname, chrom)
        #sys.stdout.flush()
    except Exception as e:
        #print(e, setname, chrom)
        #print('attempting to delete in second', setname, chrom)
        #sys.stdout.flush()
        if A is not None:
            del A
        if y is not None:
            del y
        if s is not None:
            del s
        #print('dense solver failed too. will retry later?', setname, chrom)
        #sys.stdout.flush()
        if final_try:
            #print('was final try', setname, chrom)
            gc.collect()
            #print('raising exception', setname, chrom)
            #sys.stdout.flush()
            raise Exception("Cannot allocate enough memory of solving RWR")
        else:
            #print('returning string', setname, chrom);sys.stdout.flush()
            #gc.collect()
            return "try_later"
    ##############
    #print('finalizing s', setname, chrom)
    #sys.stdout.flush()
    s *= alpha
    s += s.transpose()
    if y is not None:
        del y
    if A is not None:
        del A
    if m is not None:
        del m
    #del m, y, A
    return s


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
    #trim_value = d['value'].quantile(1-trim)
    vals = d['value'].tolist()
    vals.sort()
    vals.reverse()
    trim_value = vals[round(trim * len(vals)) - 1]
    #print(trim_value)
    remaining = d[d['value'] <= trim_value]['value']
    mu = np.mean(remaining)
    sd = np.std(remaining)
    #print(trim_value, remaining.shape[0], mu, sd)
    d.loc[:,'value'] = (d['value'] - mu) / sd
    return d

def determine_proc_share(indir, chrom_lens, n_proc, rank, outdir, ignore_sets = set(), rewrite = True):
    filenames = os.listdir(indir)
    filenames = [name for name in filenames if name.endswith(".bedpe")]
    filenames.sort()
    #print(filenames)
    setnames = [re.search(r".*\.", fname).group()[:-1] for fname in filenames]
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    chrom_list.sort(key = lambda x: x[1])
    chrom_list.reverse()
    chrom_names = [i[0] for i in chrom_list]
    #chrom_names = list(chrom_lens.keys())
    #if rank == 0:
    #    print('sorted chroms', chrom_names)
    #chrom_names.sort()
    jobs = [(chrom_names[i], filenames[j], setnames[j]) for i in range(len(chrom_names)) for j in range(len(filenames))]
    #print('init jobs len', len(jobs), '.ignore sets len', len(ignore_sets))
    jobs = [job for job in jobs if (job[2], job[0]) not in ignore_sets]
    #print('jobs len after removing ignores', len(jobs), jobs[0])
    #print('example ignore', list(ignore_sets)[0] if len(ignore_sets) > 0 else "0")
    #print('init jobs len', len(jobs))
    if rewrite:
        completed_filenames = os.listdir(outdir)
        #print(completed_filenames[0])
        incompl = [name for name in completed_filenames if not name.endswith(".normalized.rwr.bedpe")]
        if rank == 0:
            pass
            #print(incompl)
        completed_setnames = [name[:-len(".normalized.rwr.bedpe")] for name in completed_filenames if name.endswith(".normalized.rwr.bedpe")]
        #print(completed_filenames[0])
        #print('filenames', len(completed_filenames))
        completed_filenames.sort()
        #completed_setnames = [re.search(r".*\..*?\.", fname).group()[:-1] for fname in completed_filenames]
        #print('completed_senames', len(completed_setnames))
        #print(completed_setnames[0])
        completed_pairs = set([(setname[:setname.rfind('.')], setname[(setname.rfind('.') + 1):]) for setname in completed_setnames])
        #print(len(completed_pairs), 'example completed', list(completed_pairs)[0] if len(completed_pairs) > 0 else "0")
        #print('jobs vs complted', len(jobs), len(completed_pairs))
        #print(completed_pairs[:2])
        #print(jobs[:2])
        jobs = [job for job in jobs if (job[2], job[0]) not in completed_pairs]
        #print('new jobs len', len(jobs))
    jobs.sort()
    indices = list(range(rank, len(jobs), n_proc))
    #random.seed(4)
    #random.shuffle(jobs)
    #random.seed()
    #if rank in [1,2,3, 19]:
    #    print(rank, [(i[0], i[2]) for i in jobs[:5]])
    proc_jobs = [jobs[i] for i in indices]
    #proc_jobs = deque(proc_jobs)
    #proc_jobs.rotate(rank)
    #proc_jobs = list(proc_jobs)
    random.shuffle(proc_jobs)
    #print('jobs for rank', rank, len(proc_jobs), 'from total', len(jobs));print(n_proc, rank);
    return proc_jobs
    
def ammend_ignore_list(logfile, ignore_filename):
    ignore = []
    if logfile:
        with open(logfile, 'r') as lfile:
            lines = lfile.readlines()
        lines = [line.split() for line in lines]
        started = set()
        completed = set()
        for line in lines:
            if len(line) < 5:
                continue
            if line[0] == "solved":
                #remove .bedpe suffix from setname
                completed.add((line[-2], line[-1][:-6]))
            else:
                started.add((line[-2], line[-1]))
        ignore = started.difference(completed)
        with open(ignore_filename, 'a') as ofile:
            for chrom, setname in ignore:
            	ofile.write(f'{setname} {chrom}\n')
        
def get_ignore_list(filename):
    with open(filename) as ifile:
        lines = ifile.readlines()
    lines = [line.split() for line in lines]
    ignores = set()
    for line in lines:
        ignores.add((line[0], line[1]))
    return ignores

def get_nth_diag_indices(mat, offset):
    rows, cols_orig = np.diag_indices_from(mat)
    cols = cols_orig.copy()
    if offset > 0:
        cols += offset
        rows = rows[:-offset]
        cols = cols[:-offset]
    return rows, cols

def normalize_along_diagonal_from_numpy(d, chrom, max_bin_distance, output_filename, binsize, remove_bins, trim = 0.01):
    #df_all = pd.DataFrame()
    if os.path.exists(output_filename):
        os.remove(output_filename)
    with open(output_filename, "a") as f:
        for offset in range(1, max_bin_distance + 1):
            r, c = get_nth_diag_indices(d, offset)
            vals_orig = d[r,c].tolist()

            vals = vals_orig.copy()
            vals.sort()
            vals.reverse()
            trim_value = vals[(round(trim * len(vals)) - 1)]
            trim_index = round(trim * len(vals)) - 1
            remaining = vals[(trim_index):]
            mu = np.mean(remaining)
            sd = np.std(remaining)
            vals_orig = (np.array(vals_orig) - mu) / sd
            df = pd.DataFrame({'x1': r, 'y1': c, 'v': vals_orig})
            df['x1'] = ((df['x1'] + remove_bins) * binsize).astype(int)
            df['y1'] = ((df['y1'] + remove_bins) * binsize).astype(int)
            df['x2'] = (df['x1'] + binsize).astype(int)
            df['y2'] = (df['y1'] + binsize).astype(int)
            df['chr1'] = chrom
            df['chr2'] = chrom
            df = df[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'v']]
            df.to_csv(f, mode='a', header=False, sep = "\t", index = False)
        #df_all = pd.concat([df_all, df], axis = 0)
    #df_all = df_all.sort_values(['x1', 'y1'])
    #df_all.to_csv(output_filename, header = False, sep = "\t", index = False)

def get_rwr_for_all(indir, outdir = None, binsize = BIN, alpha = ALPHA, dist = DIST, chrom_lens = CHROM_DICT, 
                normalize = False, n_proc = 1, rank = 0, genome = 'mouse', filter_file = None, parallel = False, 
                                rwr_logfile = None, rwr_logfilename = None, threaded_lock = None, logger = None):
    if logger:
        logger.set_rank(rank)
    if not outdir:
        outdir = indir
        outdir = os.path.join(outdir, "rwr")
    try:
        os.makedirs(outdir)
    except:
        pass
    if rwr_logfile is None:
        rwr_logfile = open(rwr_logfilename, 'a') if rwr_logfilename else None
    #print('rwr rank', rank)
    #ignore_filename = os.path.join(outdir, 'ignore_sets.txt')
    #ammend_ignore_list(rwr_logfilename, ignore_filename)
    #ignore_sets = get_ignore_list(ignore_filename)
    processor_jobs = determine_proc_share(indir, chrom_lens, n_proc, rank, outdir)#, ignore_sets = ignore_sets)
    retry_filename = os.path.join(outdir, ('_'.join([str(rank), "retry", "instances"]) + ".txt"))
    attempt_counter = 0
    attempts_allowed = 10
    #print(rank, [(name[0], name[2]) for name in processor_jobs])
    #sys.stdout.flush()
    logger.write(f'\tprocessor {rank}: {len(processor_jobs)} jobs assigned to processor {rank}.', \
                             append_time = False, allow_all_ranks = True, verbose_level = 2)
    while len(processor_jobs) > 0 and attempt_counter < attempts_allowed:
        #print(rank, 'has len', len(processor_jobs), attempt_counter)
        #sys.stdout.flush()
        for chrom, filename, setname in processor_jobs:
            logger.flush()
            gc.collect()
            #print('in for loop', rank)
            #sys.stdout.flush()
            msg = f'rank {rank}: running {chrom} {setname}\n'
            if rwr_logfile and parallel:
                rwr_logfile.Write_shared(msg.encode('utf-8'))
                rwr_logfile.Sync()
            elif rwr_logfile:
                if threaded_lock:
                    threaded_lock.acquire()
                rwr_logfile.write(msg)
                rwr_logfile.flush()
                if threaded_lock:
                    threaded_lock.release()
            #sys.stdout.flush()
            filepath = os.path.join(indir, filename)
            final_try = False if attempt_counter < attempts_allowed else True
            #print('calling rwr for set')
            #sys.stdout.flush()
            d = get_rwr_en(filepath, binsize = binsize, distance = dist, chrom = chrom, chrom_len = chrom_lens[chrom], 
                alpha = alpha, logfile = rwr_logfile, parallel = parallel, threaded_lock = threaded_lock)
            last_bin = chrom_lens[chrom] // binsize
            if isinstance(d, str) and d == "try_later":
                logger.write(f'\tprocessor {rank}: processing of {setname} {chrom} failed, likely due to lack of memory. ' +\
                              f'Adding it to {("_".join([str(rank), "retry", "instances"]) + ".txt")} file to try again later', \
                             append_time = False, allow_all_ranks = True, verbose_level = 3)
                #print ("attempting to add to the remainder")
                #sys.stdout.flush()
                with open(retry_filename, 'a') as ofile:
                    ofile.write("\t".join([chrom, filename, setname]) + "\n")
                #print('written to file',chrom, setname, "from", rank)
                #sys.stdout.flush()
                continue
            output_filename = os.path.join(outdir, ".".join([setname, chrom, "rwr", "npy"]))
            #df.sort_values(['x1', 'y1'], inplace = True)
            #df.to_csv(output_filename, sep = "\t", header = None, index = False)
            np.save(output_filename, d)
            logger.write(f'\tprocessor {rank}: RWR solved for {setname} {chrom}. Going to normalize if requested', \
                              append_time = False, allow_all_ranks = True, verbose_level = 3)
            if normalize:
                max_bin_distance = int(dist // binsize)
                remove_bins = int(50000 // binsize)
                #print(remove_bins, type(remove_bins))
                d = d[remove_bins:-remove_bins, remove_bins:-remove_bins]
                #df = df[(df['x1'] >= 50000) & (df['y1'] <= last_bin * binsize - 50000)]
                output_filename = os.path.join(outdir, ".".join([setname, chrom, "normalized", "rwr", "bedpe"]))
                normalize_along_diagonal_from_numpy(d, chrom, max_bin_distance, output_filename, binsize, remove_bins)
                del d
                gc.collect()
                df = pd.read_csv(output_filename, header = None, sep = "\t")
                df.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'v']
                #df = df.groupby(df['y1'] - df['x1'], as_index = False).apply(normalize_along_diagonal).reset_index(drop = True)
                df.sort_values(['x1', 'y1'], inplace = True)
                df.to_csv(output_filename, sep = "\t", header = None, index = False)
                del df
        if os.path.exists(retry_filename):
            logger.write(f'\tprocessor {rank}: Attempting to re-run failed jobs. Attempt: {attempt_counter + 1}/{attempts_allowed}', \
                              append_time = False, allow_all_ranks = True, verbose_level = 3)
            #print("rank", rank, ": attempting to rerun failed jobs. Attempt #", attempt_counter + 1)
            #sys.stdout.flush()
            with open(retry_filename, 'r') as infile:
                jobs = infile.readlines()
            jobs = [line.split() for line in jobs]
            processor_jobs = jobs
            os.remove(retry_filename)
            attempt_counter += 1
        else:
            processor_jobs = []
            logger.write(f'\tprocessor {rank}: finished processing my share (RWR)', \
                              append_time = False, allow_all_ranks = True, verbose_level = 3)
            #print("rank", rank, ": no remaining jobs or parser failed")
            #sys.stdout.flush()
    if attempt_counter == attempts_allowed:
        logger.write(f'\tprocessor {rank}: Failed to finish assigned jobs after {attempts_allowed} attempts.')
        #return df
    if rwr_logfile:
        rwr_logfile.close()

if __name__ == "__main__":
    pass
    #get_rwr_for_all(INDIR, normalize = True)
