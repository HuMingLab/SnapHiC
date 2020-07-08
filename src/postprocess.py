import pandas as pd
import numpy as np
import scipy as sp
import networkx as nx
import os
import glob
from scipy import ndimage
from scipy.spatial import distance
from collections import defaultdict
import subprocess

def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_names = list(chrom_lens.keys())
    chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms

def create_local_filter(base, difference):
    base[difference:-difference, difference:-difference] = 1
    return base

def create_donut_filter(base, center):
    base[center, :] = 0
    base[:, center] = 0
    base[center:,:-center] = 0
    return base
    
def create_lower_left_filter(base, center):
    base[:-center, :center] = 1
    base[:-center,:] = 0
    base[:, center:] = 0
    return base 
    
def create_filters(gap_large, gap_small):
    width = gap_large * 2 + 1
    base = np.zeros((width, width))
    difference = gap_large - gap_small
    center = gap_large
    local_filter = create_local_filter(base.copy(), difference)
    local_inverse = abs(local_filter.copy() - 1)
    lower_left_filter = create_lower_left_filter(local_inverse.copy(), center)
    donut_filter = create_donut_filter(local_inverse.copy(), center)
    horizontal_filter = local_inverse[center, :].copy().reshape((1, width))
    vertical_filter = local_inverse[:,center].copy().reshape((width, 1))
    filter_dict = {'local': local_filter, 'lower_left': lower_left_filter, \
                   'donut':donut_filter, 'horizontal': horizontal_filter, \
                   'vertical': vertical_filter}
    return filter_dict

def determine_dense_matrix_size(dist, binsize, max_mem, num_filters):
    max_mem_floats = max_mem * 1e9 
    max_mem_floats /= 8
    square_cells = max_mem_floats // num_filters
    max_dist_bins = dist//binsize
    mat_size = min(int(np.floor(np.sqrt(square_cells))), int(max_dist_bins * 10))
    if mat_size < max_dist_bins:
        raise "Specified " + str(max_mem) + "GB is not enough for constructing dense matrix with distance " + str(dist) + "."
    return mat_size


def convert_sparse_dataframe_to_dense_matrix(d, mat_size, dist, binsize, upper_limit, chrom_size):
    d.loc[:,'i'] = (d.loc[:,'x1'] // binsize).astype(int)
    d.loc[:,'j'] = (d.loc[:,'y1'] // binsize).astype(int)
    max_distance_bin = dist // binsize
    chrom_bins = int(chrom_size // binsize)
    for i in range(0, chrom_bins + 1, int(mat_size - max_distance_bin)):
        matrix_upper_bound = max(0, i - upper_limit)
        matrix_lower_bound = min(i + mat_size + upper_limit, chrom_bins + 1)
        #print("sending mat from", matrix_upper_bound, matrix_lower_bound)
        portion = d[(d['i'] >= matrix_upper_bound) & \
                    (d['j'] < matrix_lower_bound)]
        dense_mat = sp.sparse.csr_matrix((portion['outlier_count'], \
                                        ((portion['i'] - matrix_upper_bound), \
                                         (portion['j'] - matrix_upper_bound))), \
                                       shape = (matrix_lower_bound - matrix_upper_bound, \
                                                matrix_lower_bound - matrix_upper_bound))
        dense_mat = np.array(dense_mat.todense())
        '''
        if matrix_upper_bound == 0:
            pad_size = abs(i - upper_limit)
            dense_mat = np.pad(dense_mat, ((pad_size, 0), (pad_size, 0), (0, 0)), mode = 'constant')
        if matrix_lower_bound == chrom_bins:
            pad_size = upper_limit #- 1
            mat_3d = np.pad(mat_3d, ((0, pad_size), (0, pad_size), (0, 0)), mode = 'constant')
        '''
        yield dense_mat

def filtered_mean(x):
    count = np.sum(x>-1)
    if count == 0:
        return 0
    else:
        #print(np.average(x, weights = (x>-1)))
        return np.average(x, weights = (x>-1))

def apply_all_filters(mat, footprints, candidates, start_index, max_distance_bin, gap_large):
    #print("processing from ", start_index)
    #print(mat.shape)
    #print("sub between:", start_index + gap_large, max(0, start_index) + mat.shape[0])
    sub_candidates = candidates[(candidates['i'] >= start_index + gap_large) & \
                                (candidates['j'] < max(0, start_index) + mat.shape[0])].copy()
    #print(sub_candidates.shape)
    if sub_candidates.shape[0] > 0:
        sub_candidates.loc[:,'temp_i'] = (sub_candidates['i'] - max(0, start_index)).astype(int)
        sub_candidates.loc[:,'temp_j'] = (sub_candidates['j'] - max(0, start_index)).astype(int)
        for name, footprint in footprints.items():
            #print(name)
            filter_out = sp.ndimage.generic_filter(mat, function = filtered_mean, \
                                                   footprint = footprint, mode = 'constant', cval = -1)
            sub_candidates.loc[:,name] = filter_out[sub_candidates['temp_i'], sub_candidates['temp_j']]
        return sub_candidates
    else:
        empty_dict = {i:[] for i in list(sub_candidates.columns) + list(footprints.keys())}
    return pd.DataFrame(empty_dict)

def find_candidates(indir, outdir, proc_chroms, chrom_lens, fdr_thresh, gap_large, gap_small, candidate_lower_thresh, \
                    candidate_upper_thresh, binsize, dist, max_mem, num_cells, case_to_control_diff_threshold, \
                    circle_threshold_mult, donut_threshold_mult, lower_left_threshold_mult, \
                    horizontal_threshold_mult, vertical_threshold_mult, outlier_threshold_mult, filter_file):
    footprints = create_filters(gap_large, gap_small)
    matrix_max_size = determine_dense_matrix_size(dist, binsize, max_mem, len(footprints))
    for chrom in proc_chroms:
        infile = os.path.join(indir, ".".join(["significances", chrom, "bedpe"]))
        d = pd.read_csv(infile, sep = "\t")
        candidates = d[(d['y1'] - d['x1'] <= candidate_upper_thresh) & \
                   (d['y1'] - d['x1'] >= candidate_lower_thresh) & \
                   (d['case_avg'] > 0) & \
                   (d['case_avg'] - d['control_avg'] > case_to_control_diff_threshold) & \
                   (d['fdr_dist'] < fdr_thresh) &\
                   (d['outlier_count'] > outlier_threshold_mult*num_cells)]
        #print("candidates")
        #print(candidates.shape)
        #print(fdr_thresh)
        #print(num_cells)
        candidates.loc[:,'i'] = candidates['x1'] // binsize
        candidates.loc[:,'j'] = candidates['y1'] // binsize

        #print(candidates.describe())
        submatrices = convert_sparse_dataframe_to_dense_matrix(d, matrix_max_size, dist, \
                                                               binsize, gap_large, chrom_lens[chrom])
        max_distance_bin = dist // binsize
        results = []
        for i, submatrix in enumerate(submatrices):
            start_index = i * (matrix_max_size-max_distance_bin) - gap_large
            submat_result = apply_all_filters(submatrix, footprints, candidates, start_index, \
                                              max_distance_bin, gap_large)
            results.append(submat_result)
        results = pd.concat(results, axis = 0, sort = False)
        try:
            results.drop(['i', 'j', 'temp_i', 'temp_j'], axis =1, inplace = True)
            results['circle'] = results[['lower_left', 'donut', 'horizontal', 'vertical']].sum(axis = 1)
        except:
            #if there are no candidates found for any of the iterations, there won't be any temp_i, temp_j
            pass
        columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'outlier_count', 'pvalue', \
                   'fdr_chrom', 'fdr_dist', 'case_avg', 'control_avg'] + list(footprints.keys()) + \
                  ['circle']
        results = results[columns] 
        results.to_csv(os.path.join(outdir, ".".join(["nofilter", chrom, "bedpe"])), sep = "\t", index = False)
        results = results[(results['outlier_count'] > results['circle'] * circle_threshold_mult) & \
                          (results['outlier_count'] > results['donut'] * donut_threshold_mult) & \
                          (results['outlier_count'] > results['lower_left'] * lower_left_threshold_mult) & \
                          (results['outlier_count'] > results['horizontal'] * horizontal_threshold_mult) & \
                          (results['outlier_count'] > results['vertical'] * vertical_threshold_mult)]
        if filter_file:
            filter_regions = pd.read_csv(filter_file, sep = "\t", header = None)
            filter_regions.rename({0:'chr', 1:'start'}, axis = 1, inplace = True)
            for side in ['x1', 'y1']:
                results = results.merge(filter_regions, left_on = ['chr1', side], \
                                     right_on = ['chr', 'start'], how = "outer", indicator = True)
                results = results[results['_merge'] == 'left_only'].drop('_merge', axis = 1)
            results = results[columns]
        results.to_csv(os.path.join(outdir, ".".join(["candidates", chrom, "bedpe"])), sep = "\t", index = False)
        
def combine_postprocessed_chroms(directory):
    headers = "\t".join(["chr1","x1","x2","chr2","y1","y2","outlier_count",\
                         "pvalue","fdr_chrom","fdr_dist","case_avg","control_avg",\
                         "local","circle","lower_left","donut","horizontal","vertical",\
                         "cluster","cluster_type","cluster_size","neg_log10_fdr","summit"])
    output_filename_temp = os.path.join(directory, "combined.postprocessed.bedpe.temp")
    output_filename = os.path.join(directory, "combined.postprocessed.bedpe")
    input_filepattern = directory + '/clustered.candidates.*.bedpe'
    proc = subprocess.Popen("awk 'FNR>1' " + input_filepattern + ' > ' + output_filename_temp, shell = True)
    proc.communicate()
    with open(output_filename, 'w') as ofile:
        ofile.write(headers + "\n")
    proc = subprocess.Popen(" ".join(["cat",output_filename_temp,">>",output_filename]), shell = True)
    proc.communicate()
    d = pd.read_csv(output_filename, sep = "\t")
    d = d[d['summit'] == 1]
    summits_filename = os.path.join(directory, "combined.postprocessed.summits.bedpe")
    d.to_csv(summits_filename, sep = "\t", index = False)

def label_propagate(dist_matrix_binary):
    temp_labels = np.argmax(dist_matrix_binary, axis = 1)
    peak_to_label = {i: temp_labels[i] for i in range(dist_matrix_binary.shape[0])}
    label_to_peaks = defaultdict(set)
    for label in range(dist_matrix_binary.shape[0]):
        label_to_peaks[label] = set(np.where(dist_matrix_binary[label,:])[0])
    prev_label_to_peaks = {}
    while prev_label_to_peaks != label_to_peaks:
        prev_label_to_peaks = label_to_peaks.copy()
        for label, peaks in prev_label_to_peaks.items():
            labels = [peak_to_label[peak] for peak in peaks]
            if len(labels) > 0:
                label_candidate = min(labels)
                if label_candidate < label:
                    label_to_peaks[label_candidate] = label_to_peaks[label_candidate].union(peaks)
                    del label_to_peaks[label]
                    for peak in peaks:
                        peak_to_label[peak] = label_candidate
    return label_to_peaks

def cluster_peaks(outdir, proc_chroms, clustering_gap, binsize, summit_gap):
    for chrom in proc_chroms:
        input_filename = os.path.join(outdir, ".".join(["candidates", chrom, "bedpe"]))
        d = pd.read_csv(input_filename, sep = "\t")
        if d.shape[0] > 0:
            #compute pairwise distances
            d.loc[:,'i'] = (d.loc[:,'x1'] // binsize).astype(int)
            d.loc[:,'j'] = (d.loc[:,'y1'] // binsize).astype(int)
            points = d[['i', 'j']].to_numpy()
            dists = distance.cdist(points, points, 'cityblock')
            np.fill_diagonal(dists, np.Inf)
            min_dists = dists.min(axis = 0)

            #separate singletons from cluster peaks
            singleton_indices = np.where(min_dists > clustering_gap)[0]  
            singletons = d.iloc[singleton_indices, :]
            singletons.reset_index(drop = True, inplace = True)
            clusters = d.drop(singleton_indices, axis = 0).reset_index(drop = True)

            #find labels for cluster peaks
            points = clusters[['i', 'j']].to_numpy()
            dists = distance.cdist(points, points, 'cityblock')
            dists = dists <= clustering_gap
            label_to_peaks = label_propagate(dists)

            #form cluster_name column
            singletons.loc[:,'cluster'] = list(singletons.index)
            singletons.loc[:,'cluster'] = 'singleton_' + singletons['cluster'].astype(str)
            singletons['cluster_type'] = 'singleton'
            singletons['cluster_size'] = 1
            singletons['neg_log10_fdr'] = -np.log10(singletons['fdr_chrom'])
            singletons['summit'] = 1

            if clusters.shape[0] > 0:
                clusters.loc[:,'cluster'] = 0

                for counter, (label, indices) in enumerate(label_to_peaks.items()):
                    clusters.loc[list(indices), 'cluster'] = "cluster_" + str(counter)
                def compute_cluster_stats(df):
                    df['cluster_size'] = df.shape[0]
                    df['neg_log10_fdr'] = np.sum(-np.log10(df['fdr_chrom']))
                    df['summit'] = 0
                    #df.loc[df['fdr_chrom'].idxmin(),'summit'] = 1
                    return df
                clusters = clusters.groupby('cluster').apply(compute_cluster_stats)
            
                #find broad and sharp peaks
                clusters = clusters.sort_values('neg_log10_fdr', axis = 0).reset_index(drop = True)
                temp = pd.DataFrame({'row': list(range(1, clusters.shape[0] + 1)), 'nlfdr': clusters['neg_log10_fdr']})
                temp.loc[:,'row'] = temp.loc[:,'row'] / max(temp['row'])
                temp.loc[:,'nlfdr'] = temp.loc[:,'nlfdr'] / max(temp['nlfdr'])
                rows = temp['row'].copy()
                temp.loc[:,'row'] = 1/np.sqrt(2) * rows + 1/np.sqrt(2) * temp['nlfdr']
                temp.loc[:,'nlfdr'] = -1/np.sqrt(2) * rows + 1/np.sqrt(2) * temp['nlfdr']
                temp['new_row'] = list(range(temp.shape[0]))
                ref_point = temp[temp['nlfdr'] == min(temp['nlfdr'])]['new_row'].iloc[0]
                ref_value = clusters.iloc[ref_point,:]['neg_log10_fdr']
                clusters.loc[clusters['neg_log10_fdr'] < ref_value, 'cluster_type'] = 'SharpPeak'
                clusters.loc[clusters['neg_log10_fdr'] >= ref_value, 'cluster_type'] = 'BroadPeak'
            
                def find_cluster_summits(df, summit_gap):
                    summits = pd.DataFrame({i:[] for i in list(df.columns)})
                    #print(summits.shape)
                    #print(df)
                    while (df.shape[0] > 0):
                        summit = pd.DataFrame([df.loc[df['fdr_chrom'].idxmin(),:]], columns = list(df.columns))
                        #print(summit)
                        #print(summit.shape)
                        #print(type(summit))
                        #print(df.shape)
                        summits = pd.concat([summits, summit], axis = 0)
                        #print('summit shape:')
                        #print(summit.shape)
                        #print('summits shape:')
                        #print(summits.shape)
                        df = df[(abs(df['x1'] - summit.iloc[0,:]['x1']) > summit_gap) & \
                                (abs(df['y1'] - summit.iloc[0,:]['y1']) > summit_gap)]
                    return summits
                
                summits = clusters.groupby('cluster').apply(find_cluster_summits, summit_gap = summit_gap)
                summits.to_csv(os.path.join(outdir, ".".join(["clustered", "candidates", chrom, "bedpe"])), \
                               sep = "\t", index = False)
            singletons.to_csv(os.path.join(outdir, ".".join(["singletons", "candidates", chrom, "bedpe"])), \
                              sep = "\t", index = False)
            
            #final = pd.concat([singletons, clusters], axis = 0, sort = False)
            #final.drop(["i", "j"], axis = 1, inplace = True)
            #final.to_csv(os.path.join(outdir, ".".join(["clustered", "candidates", chrom, "bedpe"])), sep = "\t", index = False)

def postprocess(indir, outdir, chrom_lens, fdr_thresh, gap_large, gap_small, candidate_lower_thresh, \
                    candidate_upper_thresh, binsize, dist, clustering_gap, rank, n_proc, max_mem, \
                    num_cells, case_to_control_diff_threshold, circle_threshold_mult, donut_threshold_mult, \
                    lower_left_threshold_mult, horizontal_threshold_mult, vertical_threshold_mult, \
                    outlier_threshold_mult, filter_file, summit_gap):
    try:
        os.makedirs(outdir)
    except:
        pass
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    find_candidates(indir, outdir, proc_chroms, chrom_lens, fdr_thresh, gap_large, gap_small, candidate_lower_thresh, \
                    candidate_upper_thresh, binsize, dist, max_mem, num_cells, case_to_control_diff_threshold, \
                    circle_threshold_mult, donut_threshold_mult, lower_left_threshold_mult, \
                    horizontal_threshold_mult, vertical_threshold_mult, outlier_threshold_mult, filter_file)
    cluster_peaks(outdir, proc_chroms, clustering_gap, binsize, summit_gap)
   
