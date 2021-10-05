
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
def get_rwr_en(edge_filename, binsize = BIN, distance = DIST, chrom = list(CHROM_DICT.keys())[0], \
		chrom_len = CHROM_DICT[list(CHROM_DICT.keys())[0]], \
                alpha = ALPHA, final_try = False, logfile = None, parallel = False, \
		threaded_lock = None, method = "inverse", max_iter = None, blin_partitions = 100,
                rwr_window_size = 200, rwr_step_size = 100):
    gc.collect()
    setname = edge_filename[(edge_filename.rfind('/') + 1):]
    #print('computing rwr for', setname, chrom)
    #sys.stdout.flush()
    edgelist = pd.read_csv(edge_filename, sep = "\t", header = None, names = ["chr1", "x1", "x2", "chr2", "y1", "y2"])
    edgelist = edgelist[(edgelist['chr1'] == chrom) & (edgelist['chr1'] == edgelist['chr2'])]
    edgelist = bin_matrix(edgelist, binsize)
    NUM = int(np.ceil(chrom_len / binsize))
    #print('NUM', NUM)
    
    #print('deleting edgelist', setname, chrom)
    #sys.stdout.flush()
    #print('getting stoc matrix', setname, chrom)
    #sys.stdout.flush()
    ##g = get_stochastic_matrix_from_edgelist(edges)
    #print('matrix returned', setname, chrom)
    #sys.stdout.flush()
    ##gc.collect()
    #print('solving rwr', setname, chrom)
    #sys.stdout.flush()
    msg = f'solved equation for {chrom} {setname}\n'
    if method == "sliding_window":
        r = sp.sparse.coo_matrix(((0,), ((0,), (0,))), shape = (NUM, NUM)).todense()
        #r = pd.DataFrame({'i':r.row, 'j': r.col, 'v': r.data})
        print("NUM", NUM)
        for window_start in range(0, NUM + rwr_window_size, rwr_step_size):
            window_end = window_start + rwr_window_size
            window_edges = edgelist[(edgelist['x1'] >= window_start) & (edgelist['y1'] < window_end)]
            defaults = pd.DataFrame({'x1':list(range(window_start, min(NUM-1, window_end - 1))), 'y1': list(range(window_start+1, min(NUM, window_end)))})
            if defaults.shape[0] == 0:
                print("skipping", window_edges.head())
                continue
            edges = pd.concat([defaults[['x1', 'y1']], window_edges[['x1', 'y1']]], axis = 0)
            edges.loc[:,'weight'] = 1
            edges.to_csv("temp2.csv")
            g = get_stochastic_matrix_from_edgelist(edges)
            partial_r = solve_rwr_inverse(g, alpha, final_try, setname, chrom)
            if isinstance(r, str):
                break
            partial_r = sp.sparse.coo_matrix(partial_r)
            partial_r = pd.DataFrame({'i':partial_r.row, 'j': partial_r.col, 'v': partial_r.data})
            partial_r = partial_r[(partial_r['i'] + partial_r['j'] > rwr_step_size) & (partial_r['i'] + partial_r['j'] < rwr_window_size + rwr_step_size)]
            #partial_r = partial_r[(partial_r['i'] + partial_r['j'] < rwr_window_size + rwr_step_size)]
            partial_r = partial_r[(partial_r['j'] - partial_r['i'] < rwr_window_size)]
            print("this is the partial output:", window_start)
            print(partial_r.describe())
            print(partial_r.head())
            partial_r['i'] += window_start
            partial_r['j'] += window_start
            #old_r = pd.merge(r, partial_r, on = ['i', 'j'], how = "outer", suffixes = ["", "_2"], indicator = True)
            #old_r = old_r.loc[old_r['_merge'] == "left_only", ['i', 'j', 'v']]
            #print("oldies:", old_r.shape)
            #r = pd.concat([old_r, partial_r], axis = 0)
            #print("newsies:", r.shape)
            #r['v'] = r.v.where((temp_r['_merge'] in ["both", "right_only"]), temp_r['v_y'], temp_r['v_x']) 
            r[partial_r['i'], partial_r['j']] = 0
            partial_r = sp.sparse.coo_matrix((partial_r['v'], (partial_r['i'], partial_r['j'])), shape = (NUM, NUM))
            r += partial_r 
    else:
        edges = pd.DataFrame({'x1':list(range(0, NUM-1)), 'y1':list(range(1, NUM))})
        edges = pd.concat([edges, edgelist[['x1', 'y1']]], axis = 0)
        edges.loc[:,'weight'] = 1
        if method == "inverse":
            g = get_stochastic_matrix_from_edgelist(edges)
            r = solve_rwr_inverse(g, alpha, final_try, setname, chrom)
        elif method == "iterative":
            g = get_stochastic_matrix_from_edgelist(edges)
            r = solve_rwr_iterative(g, alpha, final_try, setname, chrom, max_iter = max_iter)
        elif method == "nblin":
            g = get_laplacian_matrix_from_edgelist(edges)
            r = solve_rwr_nblin(g, alpha, final_try, setname, chrom)
        elif method == "blin":
            g = get_laplacian_matrix_from_edgelist(edges)
            r = solve_rwr_blin(g, alpha, final_try, setname, chrom, num_parts = blin_partitions, only_q1 = False)
        elif method == "q1":
            #g = get_laplacian_matrix_from_edgelist(edges)
            g = get_stochastic_matrix_from_edgelist(edges)
            r = solve_rwr_blin(g, alpha, final_try, setname, chrom, num_parts = blin_partitions, only_q1 = True)
        else:
            raise Exception("Unrecognized RWR method")
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

def get_laplacian_matrix_from_edgelist(edgelist):
    g = nx.from_pandas_edgelist(edgelist, source = 'x1', target = 'y1', edge_attr = ['weight'], create_using = nx.Graph())
    degrees = 1/np.array([g.degree(i) for i in g.nodes()])
    degrees = sp.sparse.diags(degrees, shape = (degrees.shape[0], degrees.shape[0]))
    print("deg", degrees.shape)
    m = sp.sparse.csc_matrix(nx.adjacency_matrix(g).astype(float))
    degrees = np.sqrt(degrees)
    m = degrees @ m @ degrees
    del degrees
    return m

def solve_rwr_inverse(stoch_matrix, alpha = ALPHA, final_try = False, setname = None, chrom = None):
    gc.collect()
    m = stoch_matrix*(1-alpha)
    m = m.transpose()
    y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")
    A = y - m
    try:
        s = None
        A = A.todense()
        y = y.todense()
        s = sp.linalg.solve(A, y)
    except Exception as e:
        if A is not None:
            del A
        if y is not None:
            del y
        if s is not None:
            del s
        if final_try:
            gc.collect()
            raise Exception("Cannot allocate enough memory of solving RWR")
        else:
            return "try_later"
    ##############
    s *= alpha
    s += s.transpose()
    if y is not None:
        del y
    if A is not None:
        del A
    if m is not None:
        del m
    return s

def solve_rwr_iterative(stoch_matrix, alpha = ALPHA, final_try = False, setname = None, chrom = None, max_iter = None):
    max_iter = max_iter if max_iter else float("inf")
    gc.collect()
    m = stoch_matrix
    y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")
    s = y
    delta = float("inf")
    A = y
    #print(m.todense()[:4,:4])
    counter = 0
    while delta > 1e-6 and counter < max_iter:
        print(delta, counter)
        counter += 1
        Aold = A.copy()
        A = (1-alpha) * m * Aold + (alpha) * y
        delta = (abs(A - Aold)).max()
    #print(delta)
    A += A.transpose()
    return A

def solve_rwr_nblin(stoch_matrix, alpha = ALPHA, final_try = False, setname = None, chrom = None):
    gc.collect()
    m = stoch_matrix
    y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")

    #svd decomposition
    U, s,  V = sp.linalg.svd(m.todense())

    #eigenvalue decompositionn
    #s, U = sp.linalg.eig(m.todense())
    #V = np.transpose(U)

    s = sp.sparse.spdiags(s, 0, m.shape[0], m.shape[0], format = "csc").todense()
    lmd = sp.linalg.inv(sp.linalg.inv(s) - (1 - alpha) * V @ U)
    A = (alpha) * (y + (1 - alpha) * U @ lmd @ V)

    #if eigenvalue decompose, might not need these lines:
    A = A / np.sum(A, axis = 0)
    A += A.transpose()
    return A

def solve_rwr_blin(stoch_matrix, alpha = ALPHA, final_try = False, setname = None, chrom = None, num_parts = None, only_q1 = False):
    num_parts = num_parts if num_parts else 100
    parts = get_partitions(stoch_matrix, num_parts)
    ordering, rev_ordering = get_orders(parts)
    stoch_matrix = stoch_matrix[:, ordering]
    q1 = construct_q1(stoch_matrix, parts, 1-alpha)
    if only_q1:
        res = q1
    else:
        w2 = construct_w2(stoch_matrix, parts)
        U, s, V = sp.linalg.svd(w2)
        s = sp.sparse.spdiags(s, 0, stoch_matrix.shape[0], stoch_matrix.shape[0], format = "csc").todense()
        lmd = np.linalg.inv(np.linalg.inv(s) - (1 - alpha) * V @ q1 @ U)
        res = (alpha)*(q1 + (1 - alpha) * q1 @ U @ lmd @ V @ q1)
    res += res.transpose()
    res = res[:, rev_ordering]
    return res

def get_partitions(adj_matrix, k):
    g = nx.from_numpy_matrix(adj_matrix.todense()) 
    import metispy 
    cost, parts = metispy.part_graph(g, k)
    ps  =  {}
    counter = 0
    for i in parts:
        if i in ps:
            continue
        else:
            ps[i] = counter
            counter += 1
    l = []
    for i in parts:
        l.append(ps[i])
    return l

def construct_q1(m, parts, alpha):
    q1 =  np.zeros(shape = m.shape)
    parts = np.array(parts)
    m = m.todense()
    for i in range(max(parts) + 1):
        nodes = list(np.where(parts == i)[0])
        start = min(nodes)
        end = max(nodes) + 1
        submat = m[start:end, start:end]
        I = np.eye(submat.shape[0])
        q1[start:end, start:end] = np.linalg.inv(I - alpha * submat)
    return q1

def construct_w2(m, parts):
    w2 =  np.zeros(shape = m.shape)
    parts = np.array(parts)
    m = m.todense()
    for i in range(max(parts) + 1):
        nodes = list(np.where(parts == i)[0])
        start = min(nodes)
        end = max(nodes) + 1
        w2[start:end, :start] = m[start:end, :start]
        w2[start:end, end:] = m[start:end, end:]
        w2[:start, start:end] = m[:start, start:end]
        w2[end:, start:end] = m[end:, start:end]
    return w2

def get_orders(parts):
    parts = np.array(parts)
    ordering = np.zeros(len(parts))
    rev_ordering = np.zeros(len(parts))
    column_counter = 0
    for part_id in range(len(set(parts))):
        group_members = np.where(parts == part_id)[0]
        for column in group_members:
            ordering[column_counter] = column
            rev_ordering[column] = column_counter
            column_counter += 1
    rev_ordering = np.array(rev_ordering, dtype = np.int)
    ordering = np.array(ordering, dtype = np.int)
    return ordering, rev_ordering

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
    setnames = [os.path.basename(fname)[:-len(".bedpe")] for fname in filenames]
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
        if False and rank == 0:
            print('jobs vs complted', len(jobs), len(completed_pairs))
            #print(completed_pairs[:2])
            print(jobs[0])
            print(list(completed_pairs)[0])
        jobs = [job for job in jobs if (job[2], job[0]) not in completed_pairs]
        if False and rank == 0:
            print('new jobs len', len(jobs))
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

def normalize_along_diagonal_from_numpy(d, chrom, max_bin_distance, output_filename, binsize, remove_bins, trim = 0.01, rwr_rank_logfile = None):
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
    if rwr_rank_logfile:
    	rwr_rank_logfile.write(f"{output_filename} completed\n")
        #df_all = pd.concat([df_all, df], axis = 0)
    #df_all = df_all.sort_values(['x1', 'y1'])
    #df_all.to_csv(output_filename, header = False, sep = "\t", index = False)

def keep_eligible_distance(d, dist, binsize):
    irange = d.shape[0]
    jrange = int(dist // binsize)
    inds = [i for i in range(irange) for j in range(i, i+jrange)]
    indps = [j for i in range(irange) for j in range(i+1, i+1+jrange)]
    keep_matrix = sp.sparse.csr_matrix(([1 for i in range(len(inds))], (inds, indps)))
    keep_matrix = keep_matrix[:,:d.shape[1]]
    #print(keep_matrix[0, 198:208])
    d = keep_matrix.multiply(sp.sparse.csr_matrix(d))
    #print(d[0, 198:208])
    return d

def get_rwr_for_all(indir, outdir = None, binsize = BIN, alpha = ALPHA, dist = DIST, chrom_lens = CHROM_DICT, 
                normalize = False, n_proc = 1, rank = 0, genome = 'mouse', filter_file = None, parallel = False, 
                                rwr_logfile = None, rwr_logfilename = None, threaded_lock = None, logger = None,
                                 keep_rwr_matrix = False, rwr_method = "inverse", blin_partitions = 100, max_iter = 100,
                                 rwr_window_size = 200, rwr_step_size = 100):
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
    rwr_rank_log_filename = os.path.join(outdir, f"{rank}_rwr_log.txt")
    rwr_rank_logfile = open(rwr_rank_log_filename, 'a')
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
                alpha = alpha, logfile = rwr_logfile, parallel = parallel, threaded_lock = threaded_lock, method = rwr_method,
                max_iter = max_iter, blin_partitions = blin_partitions, rwr_window_size = rwr_window_size, rwr_step_size = rwr_step_size)
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
            output_filename = os.path.join(outdir, ".".join([setname, chrom, "rwr"]))
            #df.sort_values(['x1', 'y1'], inplace = True)
            #df.to_csv(output_filename, sep = "\t", header = None, index = False)
            buffer_size  = 50000
            if keep_rwr_matrix:
                #print('before', d.shape)
                d = keep_eligible_distance(d, dist + buffer_size, binsize)
                #print(type(d))
                #print('after', d.shape)
                sp.sparse.save_npz(output_filename, d)
                #np.save(output_filename, d)
            logger.write(f'\tprocessor {rank}: RWR solved for {setname} {chrom}. Going to normalize if requested', \
                              append_time = False, allow_all_ranks = True, verbose_level = 3)
            if normalize:
                max_bin_distance = int((dist + buffer_size) // binsize)
                remove_bins = int(buffer_size // binsize)
                #print(remove_bins, type(remove_bins))
                #d = sp.sparse.load_npz(output_filename + ".npz")
                d = d[remove_bins:-remove_bins, remove_bins:-remove_bins]
                #df = df[(df['x1'] >= 50000) & (df['y1'] <= last_bin * binsize - 50000)]
                output_filename = os.path.join(outdir, ".".join([setname, chrom, "normalized", "rwr", "bedpe"]))
                normalize_along_diagonal_from_numpy(d, chrom, max_bin_distance, output_filename, binsize, remove_bins, rwr_rank_logfile)
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
    if rwr_logfile and not parallel:
        rwr_logfile.close()
    rwr_rank_logfile.close()

if __name__ == "__main__":
    pass
    #get_rwr_for_all(INDIR, normalize = True)
