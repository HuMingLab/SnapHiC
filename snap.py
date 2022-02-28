import argparse
import os
import sys
from src.bin_reads import bin_sets
from src.rwr import get_rwr_for_all
from src.combine_cells import combine_cells, combine_chrom_hic
from src.interaction_caller import call_interactions, combine_chrom_interactions
from src.postprocess import postprocess, combine_postprocessed_chroms, compress_and_remove_rwr_files
import glob
import pandas as pd
import multiprocessing
import src.logger
import utils.validate_rwr

def main():
    parser = create_parser()
    args = parser.parse_args()
    if args.summit_gap == -1:
        args.summit_gap = 2 * args.binsize
    #print(args.filter_file)
    #with open(args.filter_file) as ifile:
    #    lines = ifile.readlines()
    #print(lines[0])
    parallel_mode, rank, n_proc, parallel_properties = determine_parallelization_options(args.parallel, args.threaded, args.num_proc)
    if rank == 0 and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if parallel_mode == "parallel":
        parallel_properties['comm'].Barrier()
    threaded = True if parallel_mode == "threaded" else False
    logger = src.logger.Logger(f'{args.outdir}/snapHiC.log', rank = rank, verbose_threshold = args.verbose, threaded = threaded)
    #logger.set_rank(rank)
    logger.dump_args(args)
    logger.write(f'Starting operation in {parallel_mode} mode')
    #print('aall is back', rank, n_proc)
    chrom_dict = parse_chrom_lengths(args.chrom, args.chr_lens, args.genome, args.max_chrom_number)
    logger.write(f'chromosome lengths file is read. Processing {len(chrom_dict)} chromosomes.')
    #print(parallel_mode)    
    logger.flush()
    #step 1; binning
    bin_dir = os.path.join(args.outdir, "binned")
    if 'bin' in args.steps:
        logger.write('starting the binning step')
        logger.flush()
        if parallel_mode == 'nonparallel':
            bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
                     chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
                     low_cutoff = args.low_cutoff, n_proc = n_proc, rank = rank, logger = logger)
        elif parallel_mode == 'parallel':
            bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
                     chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
                     low_cutoff = args.low_cutoff, n_proc = n_proc, rank = rank, logger = logger)
            parallel_properties['comm'].Barrier()
        elif parallel_mode == 'threaded':
            params = [(args.indir, args.suffix, args.binsize, bin_dir, args.chr_columns, args.pos_columns, \
                     args.low_cutoff, n_proc, i, logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(bin_sets, params)
        logger.write("binning completed")
        logger.flush()
        #print("binned")   
        
    #step 2; RWR and normalization
    rwr_dir = os.path.join(args.outdir, "rwr")
    if 'rwr' in args.steps:
        rwr_logfilename = os.path.join(rwr_dir, "log.rwr.txt")
        logger.write(f'Starting RWR step. Additional logs for cells being processed will be written to: {rwr_logfilename}')
        logger.flush()
        try:
            os.makedirs(rwr_dir)
        except:
            pass
        if parallel_mode == 'nonparallel':
            rwr_logfile = open(rwr_logfilename, 'a')
            get_rwr_for_all(indir = bin_dir, outdir = rwr_dir, binsize = args.binsize, \
                            alpha = args.alpha, dist = args.dist, chrom_lens = chrom_dict, \
                            normalize = True, n_proc = n_proc, rank = rank, genome = args.genome, \
                            filter_file = None, parallel = False, rwr_logfile = rwr_logfile, \
                            rwr_logfilename = rwr_logfilename, threaded_lock = None, logger = logger, \
                            keep_rwr_matrix = args.keep_rwr_matrix, rwr_method = args.rwr_method, 
                            blin_partitions = args.blin_partitions, max_iter = args.max_iter)
            rwr_logfile.close()
        elif parallel_mode == 'parallel':
            parallel_properties['comm'].Barrier()
            from mpi4py import MPI
            amode = MPI.MODE_WRONLY|MPI.MODE_APPEND|MPI.MODE_CREATE
            rwr_logfile = MPI.File.Open(parallel_properties['comm'], rwr_logfilename, amode)
            rwr_logfile.Set_atomicity(True)
            get_rwr_for_all(indir = bin_dir, outdir = rwr_dir, binsize = args.binsize, \
                            alpha = args.alpha, dist = args.dist, chrom_lens = chrom_dict, \
                            normalize = True, n_proc = n_proc, rank = rank, genome = args.genome, \
                            filter_file = None, parallel = True, rwr_logfile = rwr_logfile, \
                            rwr_logfilename = rwr_logfilename, threaded_lock = None, logger = logger, \
                            keep_rwr_matrix = args.keep_rwr_matrix, rwr_method = args.rwr_method, 
                            blin_partitions = args.blin_partitions, max_iter = args.max_iter, mpi_comm = parallel_properties['comm'])
            #print(rank, 'waiting for other processes')
            rwr_logfile.Close()
            parallel_properties['comm'].Barrier()
        elif parallel_mode == 'threaded':
            #rwr_logfile = open(rwr_logfilename, 'a')
            rwr_logfile  = None
            #from multiprocessing import Lock
            m = multiprocessing.Manager()
            threaded_lock = m.Lock()
            params = [(bin_dir, rwr_dir, args.binsize, args.alpha, args.dist, chrom_dict, \
                      True, n_proc, i, args.genome, None, False, rwr_logfile, \
                      rwr_logfilename, threaded_lock, logger, args.keep_rwr_matrix, \
                      args.rwr_method, args.blin_partitions, args.max_iter) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(get_rwr_for_all, params)
            #rwr_logfile.close()
        logger.write('RWR computation completed for all cells')
        logger.flush()
        #print("rwr computed")
    
    #step 3; combine cells
    hic_dir = os.path.join(args.outdir, "hic")
    if 'hic' in args.steps:
        logger.write("Checking for completion of the RWR step")
        if rank == 0:
            rwr_is_complete = utils.validate_rwr.validate_before_combine(args.outdir, len(chrom_dict))
            if rwr_is_complete == "warning":
                logger.write("WARNING: rwr and input file counts match but there were unidentified files in the rwr directory. Please check the files listed in rwr/missing.txt to ensure their validity")
            elif not rwr_is_complete:
                print("RWR check failed! remove the files specified in rwr/missing.txt file and re-run the rwr step")
		
        logger.write("RWR check completed successfully!")
        logger.write('Combining cell-wise RWRs and generating input for juicertools pre command')
        logger.flush()
        if parallel_mode == 'nonparallel':
            combine_cells(indir = rwr_dir, outdir = hic_dir, outlier_threshold = args.outlier, \
                          chrom_lens = chrom_dict, rank = rank, n_proc = n_proc, logger = logger)
        elif parallel_mode == 'parallel':
            parallel_properties['comm'].Barrier()
            #print(rank, n_proc)
            #print('calling combine cells')
            combine_cells(indir = rwr_dir, outdir = hic_dir, outlier_threshold = args.outlier, \
                          chrom_lens = chrom_dict, rank = rank, n_proc = n_proc, logger = logger)
            parallel_properties['comm'].Barrier()
        elif parallel_mode == 'threaded':
            params = [(rwr_dir, hic_dir, args.outlier, chrom_dict, i, n_proc, logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(combine_cells, params)
        logger.write('Per chromosome combine step is completed. Next we will generate one file containing all chromosomes')
        if rank == 0:
            combine_chrom_hic(directory = hic_dir, no_cool = args.no_cool, no_hic = args.no_hic, \
                              genome = args.genome, chrom_sizes_filename = args.chr_lens, \
                              binsize = args.binsize, prefix = args.prefix)
        logger.write('Combine step is completed')
        logger.flush()
        #print("combined") 
    
    #step 4; call loops by local neighborhoods
    interaction_dir = os.path.join(args.outdir, "interactions")
    if 'interaction' in args.steps:
        logger.write('Starting to compute local background')
        logger.flush()
        if parallel_mode == 'nonparallel':
            call_interactions(indir = hic_dir, outdir = interaction_dir, chrom_lens = chrom_dict, \
                             binsize = args.binsize, dist = args.dist, \
                             neighborhood_limit_lower = args.local_lower_limit, \
                             neighborhood_limit_upper = args.local_upper_limit, rank = rank, \
                             n_proc = n_proc, max_mem = args.max_memory, logger = logger)
        elif parallel_mode == 'parallel':
            parallel_properties['comm'].Barrier()
            call_interactions(indir = hic_dir, outdir = interaction_dir, chrom_lens = chrom_dict, \
                             binsize = args.binsize, dist = args.dist, \
                             neighborhood_limit_lower = args.local_lower_limit, \
                             neighborhood_limit_upper = args.local_upper_limit, rank = rank, \
                             n_proc = n_proc, max_mem = args.max_memory, logger = logger)
            parallel_properties['comm'].Barrier()
        elif parallel_mode == 'threaded':
            params = [(hic_dir, interaction_dir, chrom_dict, args.binsize, args.dist, \
                       args.local_lower_limit, args.local_upper_limit, i, n_proc, \
                       args.max_memory, logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(call_interactions, params)
        if rank == 0:
            combine_chrom_interactions(directory = interaction_dir)
        logger.write('Local background computation is completed')
        logger.flush()
        #print("loops called")    
    
    #step 5; find candidates and cluster peaks
    postproc_dir = os.path.join(args.outdir, "postprocessed")
    if 'postprocess' in args.steps:
        logger.write('Postprocessing based on local background values')
        logger.flush()
        #num_cells = pd.read_csv(glob.glob(os.path.join(hic_dir, "*.normalized.combined.bedpe"))[0], sep = "\t").shape[1] - 7
        if parallel_mode == 'nonparallel':
            postprocess(indir = interaction_dir, outdir = postproc_dir, chrom_lens = chrom_dict, \
                        fdr_thresh = args.fdr_threshold, gap_large = args.postproc_gap_large, \
                        gap_small = args.postproc_gap_small, candidate_lower_thresh = args.candidate_lower_distance, \
                        candidate_upper_thresh = args.candidate_upper_distance, binsize = args.binsize, \
                        dist = args.dist, clustering_gap = args.clustering_gap, rank = rank, \
                        n_proc = n_proc, max_mem = args.max_memory, \
                        tstat_threshold = args.tstat_threshold, \
                        circle_threshold_mult = args.circle_threshold_multiplier, \
                        donut_threshold_mult = args.donut_threshold_multiplier, \
                        lower_left_threshold_mult = args.lower_left_threshold_multiplier, \
                        horizontal_threshold_mult = args.horizontal_threshold_multiplier, \
                        vertical_threshold_mult = args.vertical_threshold_multiplier, \
                        outlier_threshold_mult = args.outlier_threshold_multiplier, filter_file = args.filter_file, \
                        summit_gap = args.summit_gap, logger = logger)
        elif parallel_mode == 'parallel':
            parallel_properties['comm'].Barrier()
            postprocess(indir = interaction_dir, outdir = postproc_dir, chrom_lens = chrom_dict, \
                        fdr_thresh = args.fdr_threshold, gap_large = args.postproc_gap_large, \
                        gap_small = args.postproc_gap_small, candidate_lower_thresh = args.candidate_lower_distance, \
                        candidate_upper_thresh = args.candidate_upper_distance, binsize = args.binsize, \
                        dist = args.dist, clustering_gap = args.clustering_gap, rank = rank, \
                        n_proc = n_proc, max_mem = args.max_memory, \
                        tstat_threshold = args.tstat_threshold, \
                        circle_threshold_mult = args.circle_threshold_multiplier, \
                        donut_threshold_mult = args.donut_threshold_multiplier, \
                        lower_left_threshold_mult = args.lower_left_threshold_multiplier, \
                        horizontal_threshold_mult = args.horizontal_threshold_multiplier, \
                        vertical_threshold_mult = args.vertical_threshold_multiplier, \
                        outlier_threshold_mult = args.outlier_threshold_multiplier, filter_file = args.filter_file, \
                        summit_gap = args.summit_gap, logger = logger)
            parallel_properties['comm'].Barrier()
        elif parallel_mode == 'threaded':
            params = [(interaction_dir, postproc_dir, chrom_dict, args.fdr_threshold, args.postproc_gap_large, \
                        args.postproc_gap_small, args.candidate_lower_distance, \
                        args.candidate_upper_distance, args.binsize, args.dist, args.clustering_gap, \
                        i, n_proc, args.max_memory, args.tstat_threshold, \
                        args.circle_threshold_multiplier, args.donut_threshold_multiplier, \
                        args.lower_left_threshold_multiplier, args.horizontal_threshold_multiplier, \
                        args.vertical_threshold_multiplier, args.outlier_threshold_multiplier, \
                        args.filter_file, args.summit_gap, logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(postprocess, params)
        if rank == 0:
            combine_postprocessed_chroms(directory = postproc_dir, prefix = args.prefix)
            if args.compress_rwr:
                compress_and_remove_rwr_files(directory = postproc_dir)
        logger.write('Postprocessing step completed')
        logger.flush()
    logger.write('Exiting the program')
    logger.flush()

def parse_chrom_lengths(chrom, chrom_lens_filename, genome, max_chrom_number):
    if not max_chrom_number or max_chrom_number == -1:
        if not chrom or chrom == "None":
            chrom_count = 22 if genome.startswith('hg') else 19 if genome.startswith("mm") else None
            if not chrom_count:
                raise("Genome name is not recognized. Use --max-chrom-number")
            chrom = ['chr' + str(i) for i in range(1, chrom_count + 1)]
        else:
            chrom = [c.strip() for c in chrom.split()]
    else:
        chrom = ['chr' + str(i) for i in range(1, max_chrom_number + 1)]
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
            #print('it;s parallel')
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            n_proc = comm.Get_size()
            rank = comm.Get_rank()
            #print('n proc is', n_proc)
            #print('rank is', rank)
            mode = 'parallel'
            properties = {'comm': comm}
        elif threaded:
            import multiprocessing
            if n_proc < 1:
                raise Exception('if threaded flag is set, n should be a positive integer')
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
    parser.add_argument('-s', '--suffix', required = False, default = "", \
                        help = 'suffix of the input files')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-c', '--chr-columns', action = 'store', nargs = 2, \
                        type = int, help = 'two integer column numbers for chromosomes', required = False, default = "2 4")
    parser.add_argument('-p', '--pos-columns', action = 'store', nargs = 2, \
                        type = int, help = 'two integer column numbers for read positions', required = False, default = "3 5")
    parser.add_argument('-l', '--chr-lens', action = 'store', \
                        help = 'path to the chromosome lengths file', required = True)
    parser.add_argument('-g', '--genome', action = 'store', help = 'genome name; hgxx or mmxx', \
                        required = False, default = "mm10")
    parser.add_argument('--chrom', action = 'store', help = 'chromosome to process', \
                        required = False, default = None)
    parser.add_argument('--dist', type = int, help = 'distance from diagonal to consider', \
                        default = 2e6, required = False)
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
                        required = False, default = 0, type = int)
    parser.add_argument('--outlier', required = False, default = 1.96, type = float, \
                        help = 'zscore threshold for finding outliers.')
    parser.add_argument('--tstat-threshold', default = 3, type = float, required = False,
                        help = 'difference threshold for case mean to control mean for finding candidates')
    parser.add_argument('--local-lower-limit', default = 2, type = int, required = False, \
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
    parser.add_argument('--candidate-upper-distance', default = 2e6, type = int, required = False, \
                       help = 'upper threshold for distance between candidate peak binpairs')
    parser.add_argument('--circle-threshold-multiplier', default = 1.33, type = float, required = False, \
                        help = 'multiplier for circle filter threshold used in finding candidates')
    parser.add_argument('--donut-threshold-multiplier', default = 1.33, type = float, required = False, \
                        help = 'multiplier for donut filter threshold used in finding candidates')
    parser.add_argument('--lower-left-threshold-multiplier', default = 1.33, type = float, required = False, \
                        help = 'multiplier for lowerleft filter threshold used in finding candidates')
    parser.add_argument('--horizontal-threshold-multiplier', default = 1.2, type = float, required = False, \
                        help = 'multiplier for horizontal filter threshold used in finding candidates')
    parser.add_argument('--vertical-threshold-multiplier', default = 1.2, type = float, required = False, \
                        help = 'multiplier for vertical filter threshold used in finding candidates')
    parser.add_argument('--outlier-threshold-multiplier', default = 0.1, type = float, required = False, \
                        help = 'multiplier for number of outlier cells for threshold used in finding candidates')
    parser.add_argument('--summit-gap', default = 2e4, type = int, required = False,
                        help = 'disallowed distance between summit points of a cluster')
    parser.add_argument('--filter-file', default = None, required = True, \
                        help = "bed file of regions to be filtered. Regions should be binned")
    parser.add_argument('--clustering-gap', default = 1, type = int, required = False, \
                        help = 'number of allowed gaps between peaks in same cluster')
    parser.add_argument('--max-memory', default = 2, type = float, required = False, \
                        help = 'memory available in GB, that will be used in constructing dense matrices')
    parser.add_argument('--verbose', type = int, required = False, default = 0,
                        help = 'integer between 0 and 3 (inclusive), 0 for the least amount to print to log file') 
    parser.add_argument('--steps', nargs = "*", default = ['bin','rwr','hic','interaction','postprocess'], \
                        required = False, help = 'steps to run. Combination of bin,rwr,hic,interaction, and ' +\
                        'postprocess are accepted (no comma, space separated). Default is all steps.')
    parser.add_argument('--no-hic', action = 'store_true', default = False, \
                        help = 'if set, will skip generation of .hic file.', required = False)
    parser.add_argument('--no-cool', action = 'store_true', default = False, \
                        help = 'if set, will skip generation of .mcool file.', required = False)
    parser.add_argument('--keep-rwr-matrix', action = 'store_true', default = False, \
                        help = 'if set, will store the computed rwr matrix in entirety in npy format.', required = False)
    parser.add_argument('--max-chrom-number', action = "store", required = False, type = int, \
                        help = "biggest chromosome number to consider in genome, for example 22 for hg", default = -1)
    parser.add_argument('--prefix', action = 'store', required = False, default = None,
                        help = "a prefix that will be added to the name of the final output files")
    parser.add_argument('--rwr-method', action = "store", required = False, default = "inverse",
                        help = "method used for RWR computation: inverse, iterative, blin, nblin")
    parser.add_argument('--blin-partitions' , action = 'store', required = False, default = 100,
                        help = "number of METIS paritions if blin method for RWR is used.", type = int)
    parser.add_argument('--max-iter', action = 'store', required = False, type = int, default = 30,
                        help = 'maximum number of iterations if iterative method for RWR is used.')
    parser.add_argument('--compress-rwr', action = 'store_true', default = False, \
                        help = 'if set, will compress normalized rwr values and remove the original bedpe files.', required = False)
    return parser
    
if __name__ == "__main__":
    main()
