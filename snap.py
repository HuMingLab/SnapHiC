import argparse
import os
from src.bin_reads import bin_sets
from src.rwr import get_rwr_for_all
from src.combine_cells import combine_cells, combine_chrom_hic
from src.interaction_caller import call_interactions, combine_chrom_interactions
from src.postprocess import postprocess, combine_postprocessed_chroms

def main():
    parser = create_parser()
    args = parser.parse_args()
    parallel_mode, rank, n_proc, properties = determine_parallelization_options(args.parallel, args.threaded, args.num_proc)
    chrom_dict = parse_chrom_lengths(args.chrom, args.chr_lens, args.genome)
    
    
    #step 1; binning
    bin_dir = os.path.join(args.outdir, "binned")
    if parallel_mode == 'nonparallel':
        bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
                 chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
                 low_cutoff = args.low_cutoff, n_proc = n_proc, rank = rank)
    elif parallel_mode == 'parallel':
        bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
                 chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
                 low_cutoff = args.low_cutoff, n_proc = n_proc, rank = rank)
        properties['comm'].Barrier()
    elif parallel_mode == 'threaded':
        pass #with open pool, run multithreaded
    #print("binned")   
        
    #step 2; RWR and normalization
    rwr_dir = os.path.join(args.outdir, "rwr")
    if parallel_mode == 'nonparallel':
        get_rwr_for_all(indir = bin_dir, outdir = rwr_dir, binsize = args.binsize, \
                        alpha = args.alpha, dist = args.dist, chrom_lens = chrom_dict, \
                        normalize = False, n_proc = n_proc, rank = rank, genome = args.genome)
    elif parallel_mode == 'parallel':
        properties['comm'].Barrier()
        get_rwr_for_all(indir = bin_dir, outdir = rwr_dir, binsize = args.binsize, \
                        alpha = args.alpha, dist = args.dist, chrom_lens = chrom_dict, \
                        normalize = False, n_proc = n_proc, rank = rank, genome = args.genome)
        properties['comm'].Barrier()
    elif parallel_mode == 'threaded':
        pass #with open pool, run multithreaded
    #print("rwr computed")
    
    #step 3; combine cells
    hic_dir = os.path.join(args.outdir, "hic")
    if parallel_mode == 'nonparallel':
        combine_cells(indir = rwr_dir, outdir = hic_dir, outlier_threshold = args.outlier, \
                      chrom_lens = chrom_dict, rank = rank, n_proc = n_proc)
    elif parallel_mode == 'parallel':
        properties['comm'].Barrier()
        combine_cells(indir = rwr_dir, outdir = hic_dir, outlier_threshold = args.outlier, \
                      chrom_lens = chrom_dict, rank = rank, n_proc = n_proc)
        properties['comm'].Barrier()
    elif parallel_mode == 'threaded':
        pass #with open pool, run multithreaded
    if rank == 0:
        combine_chrom_hic(directory = hic_dir)
    #print("combined") 
    
    #step 4; call loops by local neighborhoods
    interaction_dir = os.path.join(args.outdir, "interactions")
    if parallel_mode == 'nonparallel':
        num_cells = call_interactions(indir = hic_dir, outdir = interaction_dir, chrom_lens = chrom_dict, \
                         binsize = args.binsize, dist = args.dist, \
                         neighborhood_limit_lower = args.local_lower_limit, \
                         neighborhood_limit_upper = args.local_upper_limit, rank = rank, \
                         n_proc = n_proc, max_mem = args.max_memory)
    elif parallel_mode == 'parallel':
        properties['comm'].Barrier()
        num_cells = call_interactions(indir = hic_dir, outdir = interaction_dir, chrom_lens = chrom_dict, \
                         binsize = args.binsize, dist = args.dist, \
                         neighborhood_limit_lower = args.local_lower_limit, \
                         neighborhood_limit_upper = args.local_upper_limit, rank = rank, \
                         n_proc = n_proc, max_mem = args.max_memory)
        properties['comm'].Barrier()
    elif parallel_mode == 'threaded':
        pass #with open pool, run multithreaded
    if rank == 0:
        combine_chrom_interactions(directory = interaction_dir)
    #print("loops called")    
    
    #step 5; find candidates and cluster peaks
    postproc_dir = os.path.join(args.outdir, "postprocessed")
    if parallel_mode == 'nonparallel':
        postprocess(indir = interaction_dir, outdir = postproc_dir, chrom_lens = chrom_dict, \
                    fdr_thresh = args.fdr_threshold, gap_large = args.postproc_gap_large, \
                    gap_small = args.postproc_gap_small, candidate_lower_thresh = args.candidate_lower_distance, \
                    candidate_upper_thresh = args.candidate_upper_distance, binsize = args.binsize, \
                    dist = args.dist, clustering_gap = args.clustering_gap, rank = rank, \
                    n_proc = n_proc, max_mem = args.max_memory, num_cells = num_cells)
    elif parallel_mode == 'parallel':
        properties['comm'].Barrier()
        postprocess(indir = interaction_dir, outdir = postproc_dir, chrom_lens = chrom_dict, \
                    fdr_thresh = args.fdr_threshold, gap_large = args.postproc_gap_large, \
                    gap_small = args.postproc_gap_small, candidate_lower_thresh = args.candidate_lower_distance, \
                    candidate_upper_thresh = args.candidate_upper_distance, binsize = args.binsize, \
                    dist = args.dist, clustering_gap = args.clustering_gap, rank = rank, \
                    n_proc = n_proc, max_mem = args.max_memory, num_cells = num_cells)
        properties['comm'].Barrier()
    elif parallel_mode == 'threaded':
        pass #with open pool, run multithreaded
    if rank == 0:
        combine_postprocessed_chroms(directory = postproc_dir)
    #print("postprocessed")

def parse_chrom_lengths(chrom, chrom_lens_filename, genome):
    if not chrom:
        chrom_count = 22 if genome == 'human' else 19 if genome == "mouse" else None
        if not chrom_count:
            raise("Genome name is not recognized")
        chrom = ['chr' + str(i) for i in range(1, chrom_count + 1)]
    else:
        chrom = [chrom]
    with open(chrom_lens_filename) as infile:
        lines = infile.readlines()
    chrom_lens = {line.split()[0]: int(line.split()[1]) for line in lines if line.split()[0] in chrom}
    return chrom_lens
    
    
def determine_parallelization_options(parallel, threaded, n_proc):
    if parallel and threaded:
        raise "Only one of 'parallel' or 'threaded' flags can be set. \
                If using a job scheduling system on a cluster with multiple machines, use parallel.\
                If using a single machine with multiple CPUs, use threaded"
    else:
        if parallel:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            n_proc = comm.Get_size()
            rank = comm.Get_rank()
            mode = 'parallel'
            properties = {'comm': comm}
        elif threaded:
            from multiprocessing import Pool
            n_proc = n_proc
            mode = 'threaded'
            rank = 0
            properties = {}
        else:
            mode = 'nonparallel'
            n_proc = 1
            rank = 0
            properties = {}
    return mode, rank, n_proc, properties
    
def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', action = 'store', required = True, \
                        help = 'input directory')
    parser.add_argument('-s', '--suffix', required = True, \
                        help = 'suffix of the input files')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-c', '--chr-columns', action = 'store', nargs = 2, \
                        type = int, help = 'two integer column numbers for chromosomes', required = True)
    parser.add_argument('-p', '--pos-columns', action = 'store', nargs = 2, \
                        type = int, help = 'two integer column numbers for read positions', required = True)
    parser.add_argument('-l', '--chr-lens', action = 'store', \
                        help = 'path to the chromosome lengths file', required = True)
    parser.add_argument('-g', '--genome', action = 'store', help = 'genome name; mouse/human', \
                        required = False, default = "mouse")
    parser.add_argument('--chrom', action = 'store', help = 'chromosome to process', \
                        required = False, default = None)
    parser.add_argument('--dist', type = int, help = 'distance from diagonal to consider', \
                        default = 1e6, required = False)
    parser.add_argument('--binsize', type = int, help = 'bin size used for binning the reads', \
                        required = False, default = 1e4)
    parser.add_argument('--low-cutoff', type = int, help = 'cut-off for removing short-range reads', \
                        default = 5e3, required = False)
    parser.add_argument('--alpha', type = float, help = 'restart probability of random walk', \
                        default = 0.05, required = False)
    parser.add_argument('--parallel', action = 'store_true', default = False, \
                        help = 'if set, will attempt to run in parallel mode', required = False)
    parser.add_argument('--threaded', action = 'store_true', default = False, \
                        help = 'if set, will attempt to use multiprocessing on single machine', required = False)
    parser.add_argument('-n', '--num-proc', help = 'number of processes used in threaded mode',
                        required = False, default = 1, type = int)
    parser.add_argument('--outlier', required = False, default = 0.97500211, type = float, \
                        help = 'percentage threshold for finding outliers.')
    parser.add_argument('--local-lower-limit', default = 3, type = int, required = False, \
                        help = 'number of bins around center (in each direction) to exlude from neighborhood')
    parser.add_argument('--local-upper-limit', default = 5, type = int, required = False, \
                        help = 'number of bins around center (in each direction) forming the neighborhood')
    parser.add_argument('--fdr-threshold', default = 0.1, type = float, required = False, \
                        help = 'FDR threshold used for candidate peak detection')
    parser.add_argument('--postproc-gap-large', default = 5, type = int, required = False, \
                       help = 'number of bins around peaks to consider in postprocessing')
    parser.add_argument('--postproc-gap-small', default = 2, type = int, required = False, \
                       help = 'number of bins around peaks to disregard in postprocessing')
    parser.add_argument('--candidate-lower-distance', default = 100000, type = int, required = False, \
                       help = 'lower threshold for distance between candidate peak binpairs')
    parser.add_argument('--candidate-upper-distance', default = 900000, type = int, required = False, \
                       help = 'upper threshold for distance between candidate peak binpairs')
    parser.add_argument('--clustering-gap', default = 3, type = int, required = False, \
                        help = 'number of allowed gaps between peaks in same cluster')
    parser.add_argument('--max-memory', default = 2, type = float, required = False, \
                        help = 'memory available in GB, that will be used in constructing dense matrices')
    return parser
    

if __name__ == "__main__":
    main()