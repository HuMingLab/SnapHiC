import os
import subprocess
import scipy
from scipy import stats
import pandas as pd
import glob
import numpy as np

def combine_and_reformat_chroms(indir, output_filename, chrom, outlier_threshold):
    input_filepattern = indir + '/*' + chrom + '*.bedpe'
    input_filenames = glob.glob(input_filepattern)
    #print(input_filenames)
    
    #copy all the rwr values from all cells for the chrom to one file
    output_file = open(output_filename, 'w')
    command = 'eval paste $(printf "<(cut -f7 %s) " ' + input_filepattern + ')'
    #print(command)
    subprocess.check_call(command, stdout = output_file, shell = True, executable = "/bin/bash")
    #proc.communicate()
    output_file.close()
    
    #count the number of outliers
    outlier_norm_count = scipy.stats.norm.ppf(outlier_threshold)
    d = pd.read_csv(output_filename, sep = "\t", header = None)
    mat = d.to_numpy()
    outliers = np.sum(mat > outlier_norm_count, axis = 1)
    
    #add the bins as the first columns of bedpe file, and append outliers as last column
    binpairs = pd.read_csv(input_filenames[0], sep = "\t", header = None).iloc[:,:-1]
    d = pd.concat([binpairs, d], axis = 1)
    d.loc[:, 'outliers'] = outliers
    
    #write output to file
    d.to_csv(output_filename, index = False, header = None, sep = "\t")
    
    #reformat
    num_cells = d.shape[1] - 7
    hic_format = d.iloc[:,:6]
    hic_format.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']
    hic_format.loc[:,'score'] = np.ceil(d.iloc[:, -1] / num_cells * 100).astype(int)
    hic_format.loc[:,'str1'] = 0
    hic_format.loc[:,'str2'] = 1
    hic_format.loc[:,'frag1'] = 0
    hic_format.loc[:,'frag2'] = 1
    hic_format = hic_format[['str1', 'chr1', 'x1', 'frag1','str2', 'chr2', 'y1', 'frag2', 'score']]
    hic_format.to_csv(output_filename + ".hic.input", index = False, header = None, sep = "\t")
    
def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_names = list(chrom_lens.keys())
    chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms
    
def combine_cells(indir, outdir, outlier_threshold, chrom_lens, rank, n_proc):
    try:
        os.makedirs(outdir)
    except:
        pass
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    for chrom in proc_chroms:
        output_filename = os.path.join(outdir, ".".join([chrom, "normalized", "combined", "bedpe"]))
        combine_and_reformat_chroms(indir, output_filename, chrom, outlier_threshold)
        
def combine_chrom_hic(directory):
    output_filename = os.path.join(directory, "allChr.hic.input")
    input_filepattern = directory + '/*.bedpe.hic.input'
    proc = subprocess.Popen('cat ' + input_filepattern + ' > ' + output_filename, shell = True)
    proc.communicate()
        