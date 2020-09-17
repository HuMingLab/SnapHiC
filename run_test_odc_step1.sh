#!/bin/bash
#PBS -q hotel
#PBS -N snHiC_s1
#PBS -l nodes=15:ppn=2
#PBS -l walltime=50:00:00
#PBS -l mem=1800GB
#PBS -o log.${PBS_JOBNAME}.${PBS_JOBID}.log
#PBS -e log.${PBS_JOBNAME}.${PBS_JOBID}.log.err
#PBS -m abe
 
cd ${PBS_O_WORKDIR}

module purge
module load intel/2018.1.163
module load intelmpi/2018.1.163
unset PYTHONPATH

indir="/projects/ps-renlab/abnousa/public_html/snapHiC/test/input/Ecker/ODC"
suffix="_indexed_contacts.txt.gz"
outdir="/projects/ps-renlab/abnousa/public_html/snapHiC/test/output/Ecker/ODC"
chrs="2 4"
pos="3 5"
chrlen="ext/hg19.chrom.sizes"
genome="hg19"
fdr_thresh=0.01
filter_file="ext/hg19_filter_regions.txt"
steps="bin rwr"

mpirun -np 30 /home/abnousa/anaconda3_2/bin/python ./snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --fdr-threshold $fdr_thresh --filter-file $filter_file --steps $steps --parallel

