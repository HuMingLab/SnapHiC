import os
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd
import skimage
import subprocess
from statsmodels.stats.multitest import multipletests
import gc
import sys
import h5py

def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    chrom_list.sort(key = lambda x: x[1])
    chrom_list.reverse()
    chrom_names = [i[0] for i in chrom_list]
    #chrom_names = list(chrom_lens.keys())
    #chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms

def combine_chrom_interactions(directory):
    headers = "\t".join(["chr1", "x1", "x2","chr2","y1","y2","outlier_count",\
                         "case_avg","control_avg","pvalue","fdr_dist","fdr_chrom"])
    output_filename_temp = os.path.join(directory, "combined_significances.bedpe.temp")
    output_filename = os.path.join(directory, "combined_significances.bedpe")
    input_filepattern = directory + '/significances.*.bedpe'
    proc = subprocess.Popen("awk 'FNR>1' " + input_filepattern + ' > ' + output_filename_temp, shell = True)
    proc.communicate()
    with open(output_filename, 'w') as ofile:
        ofile.write(headers + "\n")
    proc = subprocess.Popen(" ".join(["cat",output_filename_temp,">>",output_filename]), shell = True)
    proc.communicate()

def determine_dense_matrix_size(num_cells, dist, binsize, max_mem):
    max_mem_floats = max_mem * 1e9 
    max_mem_floats /= 8
    square_cells = max_mem_floats // num_cells
    mat_size = int(np.floor(np.sqrt(square_cells)) / 4)
    if mat_size < (dist // binsize):
        raise "Specified " + str(max_mem) + "GB is not enough for constructing dense matrix with distance " + str(dist) + "."
    return mat_size

def convert_sparse_dataframe_to_dense_matrix(d, mat_size, dist, binsize, upper_limit, num_cells, chrom_size, chrom_filename):
    d['i'] = (d.iloc[:,1] // binsize).astype(int)
    d['j'] = (d.iloc[:,4] // binsize).astype(int)
    #all_rows = set(range(d.shape[0]))
    max_distance_bin = dist // binsize
    chrom_bins = int(chrom_size // binsize)
    for i in range(0, chrom_bins + 1, int(mat_size - max_distance_bin)):
        matrix_upper_bound = max(0, i - upper_limit)
        matrix_lower_bound = min(i + mat_size + upper_limit, chrom_bins + 1)
        keeprows = list(np.where((d['i'] >= matrix_upper_bound) & (d['j'] < matrix_lower_bound))[0])
        d_portion = d.iloc[keeprows, 0:6].reset_index(drop = True)
        d_portion.columns = ['chr1','x1','x2','chr2','y1','y2']
        #print(d_portion, 'd_portions shape')
        #skiprows = all_rows.difference(keeprows)
        hdf_file = h5py.File(chrom_filename + '.cells.hdf', 'r')
        portion = hdf_file[list(hdf_file.keys())[0]]
        portion = portion[keeprows, :]
        hdf_file.close()
        if portion.shape[0] == 0:
            continue
        portion = pd.DataFrame(portion)
        #print('portions shape', portion.shape)
        portion = pd.concat([d_portion, portion], axis = 1)
        #print('concatted portion', portion.shape)
        #print(np.where(portion.isnull().sum() > 0))
        portion['i'] =  (portion.loc[:,'x1'] // binsize).astype(int)
        portion['j'] =  (portion.loc[:,'y1'] // binsize).astype(int)
        #portion_old = d[(d['i'] >= matrix_upper_bound) & \
        #            (d['j'] < matrix_lower_bound)]
        #if portion.shape[0] == 0:
        #    continue
        full_sparse = pd.DataFrame({'i': range(min(portion['i']), max(portion['j'])-1), \
                    'j': range(min(portion['i'])+1, max(portion['j']))})
        portion = portion.merge(full_sparse, on = ['i','j'], how = "outer")
        #print(portion.iloc[:,list(range(7)) + [11, 12]])
        #print("start", matrix_upper_bound, "end", matrix_lower_bound)
        dense_cells = []
        #print('here', portion.columns)
        #print(portion.head())
        #print(portion.dtypes[:20])
        #sys.stdout.flush()
        for cell_index in range(num_cells):
            cell_mat = sp.sparse.csr_matrix((portion.iloc[:, 6 + cell_index], \
                                            ((portion['i'] - matrix_upper_bound), \
                                             (portion['j'] - matrix_upper_bound))), \
                                           shape = (matrix_lower_bound - matrix_upper_bound, \
                                                    matrix_lower_bound - matrix_upper_bound))
            cell_mat = np.array(cell_mat.todense())
            cell_mat[np.tril_indices(cell_mat.shape[0],0)] = np.nan
            dense_cells.append(cell_mat)
        mat_3d = np.stack(dense_cells, axis = -1)
        if matrix_upper_bound == 0:
            pad_size = abs(i - upper_limit)
            mat_3d = np.pad(mat_3d, ((pad_size, 0), (pad_size, 0), (0, 0)), mode = 'constant', constant_values = np.nan)
        if matrix_lower_bound == chrom_bins + 1:
            pad_size = upper_limit #- 1
            mat_3d = np.pad(mat_3d, ((0, pad_size), (0, pad_size), (0, 0)), mode = 'constant', constant_values = np.nan)
        yield mat_3d, i

def get_nth_diag_indices(mat, offset):
    rows, cols_orig = np.diag_indices_from(mat)
    cols = cols_orig.copy()
    if offset > 0:
        cols += offset
        rows = rows[:-offset]
        cols = cols[:-offset]
    return rows, cols

def get_neighbor_counts_matrix(shape, gap_large, gap_small, max_distance):
    a = np.zeros(shape)
    big_width = gap_large*2 + 1
    small_width = gap_small*2 + 1
    area = big_width**2 - small_width**2
    for i in range(0, big_width):
        val = np.sum(list(range(big_width - i)))
        rows, cols = get_nth_diag_indices(a, i + 1)
        a[rows, cols] -= val
        #rows, cols = get_nth_diag_indices(a, a.shape[0]-i)
        rows, cols = get_nth_diag_indices(a, max_distance-i)
        a[rows, cols] -= val
    for  i in range(0, small_width):
        val = np.sum(list(range(small_width - i)))
        rows, cols = get_nth_diag_indices(a, i + 1)
        a[rows, cols] += val
        #rows, cols = get_nth_diag_indices(a, a.shape[0]-i)
        rows, cols = get_nth_diag_indices(a, max_distance-i)
        a[rows, cols] += val
    a += area
    return a

def compute_significances(mat, upper_limit, lower_limit, num_cells, start_index, max_distance_bin):
    gc.collect()
    #sliding window
    #print('in function')
    #sys.stdout.flush()
    big_neighborhood = skimage.util.view_as_windows(mat, (2*upper_limit+1,2*upper_limit+1,num_cells), step = 1)
    small_neighborhood = skimage.util.view_as_windows(mat, (2*lower_limit+1,2*lower_limit+1,num_cells), step = 1)
    
    #reshape to (matsize, matsize, numcells, num_neighbors+1)
    big_neighborhood = np.squeeze(big_neighborhood.reshape(big_neighborhood.shape[0],\
                                               big_neighborhood.shape[0], 1, -1, num_cells))
    small_neighborhood = np.squeeze(small_neighborhood.reshape(small_neighborhood.shape[0],\
                                               small_neighborhood.shape[0], 1, -1, num_cells))
    small_neighborhood = np.swapaxes(small_neighborhood, -2, -1)
    big_neighborhood = np.swapaxes(big_neighborhood, -2, -1)
    big_neighborhood_counts = np.sum(~np.isnan(big_neighborhood), axis = -1)
    small_neighborhood_counts = np.sum(~np.isnan(small_neighborhood), axis = -1)

    #sum on the last axis (sum of neighbors)
    big_neighborhood = big_neighborhood.sum(axis = -1)
    small_neighborhood = small_neighborhood.sum(axis = -1)
    
    #remove edge cases that are used only as neighbors
    trim_size = upper_limit - lower_limit
    small_neighborhood = small_neighborhood[trim_size:-trim_size, trim_size:-trim_size]
    small_neighborhood_counts = small_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    mat = mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    
    #local_neighborhood
    local_neighborhood = big_neighborhood - small_neighborhood
    local_neighborhood_counts = big_neighborhood_counts - small_neighborhood_counts
    #local_neighbors_count = (((upper_limit*2+1) ** 2 - (lower_limit*2+1) ** 2) - (upper_limit*2+1) + (lower_limit*2+1))/2
    del small_neighborhood, big_neighborhood, big_neighborhood_counts, small_neighborhood_counts
    gc.collect()
    
    #for each cell compute the average value over the neighborhood
    ##print('local_neighborhood shape:', local_neighborhood.shape)
    ##print('neighbors_count_mat shape: ' , neighbor_counts_matrix.shape)
    ##neighbor_counts = np.repeat(neighbor_counts_matrix[:local_neighborhood.shape[0], :local_neighborhood.shape[0], np.newaxis], \
    ##                            local_neighborhood.shape[2], axis=2)
    ##print('new shape', neighbor_counts.shape)
    ##print(len(np.where(neighbor_counts == 0)[0]))
    local_neighborhood /= local_neighborhood_counts
    del local_neighborhood_counts
    gc.collect()
    
    #compute averages over all cells for each point and each local neighborhood
    ##print(local_neighborhood.shape, mat.shape)
    pvals = stats.ttest_rel(mat, local_neighborhood, axis = 2).pvalue
    ##print(pvals.shape)
    local_neighborhood = np.mean(local_neighborhood, axis = -1)
    mat = np.mean(mat, axis = -1)
    
    #keep only upper triangle
    mat = np.triu(mat, 1)
    local_neighborhood = np.triu(local_neighborhood, 1)
    pvals = np.triu(pvals, 1)
    pvals = np.nan_to_num(pvals, nan = 1)
    
    #convert matrix to dataframe
    mat = sp.sparse.coo_matrix(mat)
    local_neighborhood = sp.sparse.coo_matrix(local_neighborhood)
    pvals = sp.sparse.coo_matrix(pvals)
    result_mat = pd.DataFrame({'i': mat.row, 'j': mat.col, 'case_avg': mat.data})
    result_neighb = pd.DataFrame({'i': local_neighborhood.row, \
                                  'j': local_neighborhood.col, \
                                  'control_avg': local_neighborhood.data})
    result_pval = pd.DataFrame({'i': pvals.row, 'j': pvals.col, 'pvalue': pvals.data})
    result = result_mat.merge(result_neighb, on = ['i', 'j'], how = "outer")
    result = result.merge(result_pval, on = ['i','j'], how = "outer")
    result.loc[:,'pvalue'] = result['pvalue'].fillna(0)
    result.loc[:,'i'] += (start_index)# + upper_limit)
    result.loc[:,'j'] += (start_index)# + upper_limit)
    result.loc[:,'i'] = result['i'].astype(int)
    result.loc[:,'j'] = result['j'].astype(int)
    result = result[result['j'] - result['i'] <= max_distance_bin]
    return result
    
def call_interactions(indir, outdir, chrom_lens, binsize, dist, neighborhood_limit_lower = 3, \
                      neighborhood_limit_upper = 5, rank = 0, n_proc = 1, max_mem = 2, logger = None):
    logger.set_rank(rank)
    try:
        os.makedirs(outdir)
    except:
        pass
    
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    #print(rank, proc_chroms)
    #sys.stdout.flush()
    for chrom in proc_chroms:
        logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
        #print(rank, chrom)
        #d = pd.read_csv(chrom_filename, sep = "\t", header = None, usecols = [0,1,2,3,4,5, num_cells + 6])
        ##command = "awk -F '\t' '{print NF; exit}' " + chrom_filename
        ##proc_output = subprocess.check_output(command, shell = True, executable = "/bin/bash")
        ##num_cells = int(proc_output) - 7
        chrom_filename = os.path.join(indir, ".".join([chrom, "normalized", "combined", "bedpe"]))
        with h5py.File(chrom_filename + ".cells.hdf", 'r') as ifile:
            num_cells = ifile[chrom].shape[1]
        logger.write(f'\tprocessor {rank}: detected {num_cells} cells for chromosome {chrom}', \
                             append_time = False, allow_all_ranks = True, verbose_level = 2)
        #print('num_cells', num_cells)
        #sys.stdout.flush()
        d = pd.read_csv(chrom_filename, sep = "\t", header = None)
        #num_cells = d.shape[1] - 7
        matrix_max_size = determine_dense_matrix_size(num_cells, dist, binsize, max_mem)
        #print(rank, matrix_max_size)
        submatrices = convert_sparse_dataframe_to_dense_matrix(d, matrix_max_size, \
                                                               dist, binsize, neighborhood_limit_upper, \
                                                               num_cells, chrom_lens[chrom], chrom_filename)
        max_distance_bin = dist // binsize
        results = []
        #print(matrix_max_size, neighborhood_limit_upper, neighborhood_limit_lower)
        #neighbor_counts_matrix = get_neighbor_counts_matrix((matrix_max_size + neighborhood_limit_upper * 2, \
        #                                                     matrix_max_size + neighborhood_limit_upper * 2), \
        #                                                 neighborhood_limit_upper, \
        #                                                 neighborhood_limit_lower, max_distance_bin)
        #print('num zeros_2d', len(np.where(neighbor_counts_matrix==0)[0]))
        #print('going in for')
        #sys.stdout.flush()
        for i, (submatrix, start_index) in enumerate(submatrices):
            logger.write(f'\tprocessor {rank}: computing background for batch {i} of {chrom}, start index = {start_index}', \
                              verbose_level = 3, allow_all_ranks = True, append_time = False)
            #print('iteration', i)
            #print('start_index', start_index)
            #sys.stdout.flush()
            if i > 0:
                limit =  i * (matrix_max_size - max_distance_bin) #- neighborhood_limit_upper
                #results[-1] = results[-1][results[-1]['i'] < limit]
                results[-1] = results[-1][results[-1]['i'] < start_index]
            #start_index = i * (matrix_max_size - max_distance_bin) - neighborhood_limit_upper
            #print(start_index)
            submat_result = compute_significances(submatrix, neighborhood_limit_upper, \
                                                  neighborhood_limit_lower, num_cells, start_index, \
                                                  max_distance_bin)
            #print('returned')
            results.append(submat_result)
        #print(rank, 'offtheloop')
        #print(rank, len(results))
        results = pd.concat(results, axis = 0)
        #print(rank, chrom, results.shape[0])
        min_index = 0
        max_index = results['j'].max()
        #print(max_index, min_index, results['i'].dtype, results['j'].dtype, neighborhood_limit_upper)
        results = results[(results['i'] >= min_index + neighborhood_limit_upper) & \
                          (results['j'] <= max_index - neighborhood_limit_upper)]
        #print(results.shape[0])

        
        def compute_fdr_by_dist(d):
            fdrs = multipletests(list(d['pvalue']), method = 'fdr_bh')[1]
            d.loc[:,'fdr_dist'] = fdrs
            return d
           
        results.reset_index(drop = True, inplace = True)
        results = results.groupby(results['j'] - results['i'], as_index = False).apply(compute_fdr_by_dist)
        results.loc[:,'fdr_chrom'] = multipletests(list(results['pvalue']), method = 'fdr_bh')[1]
        results.loc[:,'i'] = (results['i'] * binsize).astype(int)
        results.loc[:,'j'] = (results['j'] * binsize).astype(int)
       
        #print('finishing', d.shape) 
        d = d.iloc[:, list(range(7))]
        d.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'outlier_count']
        #print(d.head())
        #print(results.head())
        #d = d.merge(results, left_on = ['x1', 'y1'], right_on = ['i', 'j'], how = "outer")
        d = d.merge(results, left_on = ['x1', 'y1'], right_on = ['i', 'j'])
        #print(d.shape)
        d.drop(['i', 'j'], axis =1, inplace = True)
        logger.write(f'\tprocessor {rank}: computation for {chrom} completed. writing to file.', \
                             append_time = False, allow_all_ranks = True, verbose_level = 2)
        d.to_csv(os.path.join(outdir, ".".join(["significances", chrom, "bedpe"])), sep = "\t", index = False)   
                        
                        
