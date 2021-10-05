import os
import sys
import glob

def main(args):
    outdir = args[1]
    input_filenames = set(glob.glob(f'{outdir}/rwr/*.normalized.rwr.bedpe'))
    log_files = glob.glob(f'{outdir}/rwr/[0-9]*_rwr_log.txt')
    completed_files = []
    for logfile in log_files:
        with open(logfile) as lfile:
            llines = lfile.readlines()
            rank_files = [line.split(' ')[0] for line in llines if len(line.split(' ')) == 2]
        completed_files += rank_files
    completed_files = set(completed_files)
    incomplete_files = input_filenames.difference(completed_files)
    print(f'{len(incomplete_files)} files have been detected to have faulty rwr files. Check the missig.txt file for a list of their names.')
    with open(f'{outdir}/rwr/missing.txt', 'w') as ofile:
        ofile.write('\n'.join(list(incomplete_files)))

args = sys.argv
main(args)
