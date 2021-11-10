import os
import sys
import glob

def find_incomplete_files(outdir):
    input_filenames = glob.glob(f'{outdir}/rwr/*.normalized.rwr.bedpe')
    input_filenames = set([fname.replace("//", "/") for fname in input_filenames])
    log_files = glob.glob(f'{outdir}/rwr/[0-9]*_rwr_log.txt')
    completed_files = []
    for logfile in log_files:
        with open(logfile) as lfile:
            llines = lfile.readlines()
            rank_files = [line.split(' ')[0] for line in llines if len(line.split(' ')) == 2]
        completed_files += rank_files
    completed_files = [cfile.replace("//","/") for cfile in completed_files]
    completed_files = set(completed_files)
    incomplete_files = input_filenames.difference(completed_files)
    print(f'{len(incomplete_files)} files have been detected to have faulty rwr files. Check the missig.txt file for a list of their names.')
    with open(f'{outdir}/rwr/missing.txt', 'w') as ofile:
        ofile.write('\n'.join(list(incomplete_files)))
    return len(completed_files), len(incomplete_files)

def validate_before_combine(outdir, chroms_count):
    completes_count, incompletes_count = find_incomplete_files(outdir)
    all_bedpes = len(glob.glob(f'{outdir}/rwr/*.normalized.rwr.bedpe'))
    input_count = len(glob.glob(f'{outdir}/binned/*.bedpe')) 
    expected_count = input_count * chroms_count
    if completes_count >= expected_count: 
        if incompletes_count != 0:
            return "warning"
        return True
    else:
        raise Exception(f"RWR is not completed correctly! Please remove the files specified in rwr/missing.txt and re-run the rwr step. expected={expected_count}, completed={completes_count}!")
    
    

if __name__ == "__main__":
    args = sys.argv
    outdir = args[1]
    find_incomplete_files(outdir)
