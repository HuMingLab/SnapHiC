import os
import subprocess
import scipy
from scipy import stats
import pandas as pd
import glob
import numpy as np
import re
import sys
import h5py

def combine_and_reformat_chroms(indir, output_filename, chrom, outlier_threshold, logger, rank):
    input_filepattern = indir + '/*' + chrom + '.normalized.rwr.bedpe'
    input_filenames = glob.glob(input_filepattern)
    #print(input_filenames[0])
    #print(os.path.basename(input_filenames[0]))
    
    tempfile1 = output_filename + '.tmp1'
    #tempfile2 = output_filename + '.tmp2'
    outlier_count_tempfile = output_filename + '.tmp2'
    
    #print('opening hdfs', chrom)
    #sys.stdout.flush()
    hdf_file = h5py.File(output_filename + ".cells.hdf", 'a')
    cells_data = hdf_file.create_dataset(chrom, chunks = (100000, len(input_filenames)), shape = (0, len(input_filenames)), \
                                               maxshape = (None, len(input_filenames)))        

    batch_size = 50
    for i, batch_start in enumerate(range(0, len(input_filenames), batch_size)):
        logger.write(f'\tprocessor {rank}: combining batch {i} for {chrom}', verbose_level = 3, \
                         append_time = False, allow_all_ranks = True)
        #print("batch #", i, chrom)
        #sys.stdout.flush()
        if os.path.exists(tempfile1):
            os.remove(tempfile1)
        output_file = open(tempfile1, 'w');
        #print("fope");sys.stdout.flush();
        batch_end = min(batch_start + batch_size, len(input_filenames))
        batch_files = input_filenames[batch_start:batch_end];
        #print("fope2");sys.stdout.flush();
        command = 'eval paste $(for i in ${fnames=' + ' '.join(batch_files) + '}; do echo -n "<(cut -f7 $i) "; done)';
        #print("fope3");sys.stdout.flush();
        subprocess.check_call(command, stdout = output_file, shell = True, executable = "/bin/bash");
        #print("fope4");sys.stdout.flush();
        output_file.close();
        #print("fope5");sys.stdout.flush();
        #print("command completed", chrom)
        #sys.stdout.flush()
        vals = pd.read_csv(tempfile1, sep = "\t", header = None).to_numpy()
            
        if i == 0:
            cells_data.resize((vals.shape[0]), axis = 0)
            outliers = np.sum(vals > outlier_threshold, axis = 1).astype(int)
        else:
            outliers += np.sum(vals > outlier_threshold, axis = 1).astype(int)
        #print("vals are ready", chrom)
        #sys.stdout.flush()
        cells_data[:,batch_start:batch_end] = vals
        #print("vals transfered", chrom)
        #sys.stdout.flush();print("loopend");sys.stdout.flush();
    hdf_file.close()
    #print("hdf temp closed", chrom)
    #sys.stdout.flush()
     
    if os.path.exists(tempfile1):
        os.remove(tempfile1)
    #outliers_df = pd.DataFrame({'outliers' : outliers})
    #outliers_df.to_csv(outlier_count_tempfile, index = False, header = None, sep = "\t")
    
    #add the bins as the first columns of bedpe file, and append outliers as last column
    #binpairs = pd.read_csv(input_filenames[0], sep = "\t", header = None, usecols = [0,1,2,3,4,5]) #.iloc[:,:-1]
    #command = 'cut -f1-6 ' + input_filenames[0] + ' | paste - ' + tempfile1 + ' ' + outlier_count_tempfile
    #output_file = open(output_filename, 'w')
    #subprocess.check_call(command, stdout = output_file, shell = True, executable = "/bin/bash")
    #output_file.close()

    #d = pd.concat([binpairs, d], axis = 1)
    #print(chrom, len(outliers), d.shape)
    #print(d.columns)
    #d['outliers'] = outliers
    
    #write output to file
    #d.to_csv(output_filename, index = False,  header = None, sep = "\t")
    
    #reformat
    #get number of cells
    #command = "awk -F '\t' '{print NF; exit}' " + output_filename
    #proc_output = subprocess.check_output(command, shell = True, executable = "/bin/bash") 
    #num_cells = int(proc_output) - 7
    num_cells = len(input_filenames)
    #print("num cells:", num_cells)
    #num_cells = d.shape[1] - 7

    scores = np.ceil(outliers / num_cells * 100).astype(int)
    
    binpairs = pd.read_csv(input_filenames[0], sep = "\t", header = None, usecols = [0,1,2,3,4,5])
    binpairs.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']
    binpairs['outliers'] = outliers
    binpairs.to_csv(output_filename, sep = "\t", index = False, header = None)
    
    hic_format = binpairs.drop('outliers', axis = 1)
    hic_format.loc[:,'score'] = scores
    hic_format.loc[:,'str1'] = 0
    hic_format.loc[:,'str2'] = 1
    hic_format.loc[:,'frag1'] = 0
    hic_format.loc[:,'frag2'] = 1
    hic_format = hic_format[['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score']]
    hic_format.to_csv(output_filename + ".hic.input", index = False, header = None, sep = "\t")
    
def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    chrom_list.sort(key = lambda x: x[1])
    chrom_list.reverse()
    chrom_names = [i[0] for i in chrom_list]
    #if rank == 0:
    #    print(chrom_names)
    #chrom_names = list(chrom_lens.keys())
    #chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms
    
def combine_cells(indir, outdir, outlier_threshold, chrom_lens, rank, n_proc, logger):
    logger.set_rank(rank)
    #print('in combine cells')
    try:
        os.makedirs(outdir)
    except:
        pass
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    #print(rank, proc_chroms)
    for chrom in proc_chroms:
        logger.write(f'\tprocessor {rank}: combining chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
        #print(rank, 'processing', chrom)
        output_filename = os.path.join(outdir, ".".join([chrom, "normalized", "combined", "bedpe"]))
        combine_and_reformat_chroms(indir, output_filename, chrom, outlier_threshold, logger, rank)
        logger.write(f'\tprocessor {rank}: chromosome {chrom} is combined', verbose_level = 1, allow_all_ranks = True)
        
def combine_chrom_hic(directory):
    output_filename = os.path.join(directory, "allChr.hic.input")
    input_filepattern = directory + '/*.bedpe.hic.input'
    proc = subprocess.Popen('cat ' + input_filepattern + ' > ' + output_filename, shell = True)
    proc.communicate()
        
