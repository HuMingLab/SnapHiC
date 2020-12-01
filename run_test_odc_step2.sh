#!/bin/bash
#PBS -q hotel
#PBS -N snHiC_s2
#PBS -l nodes=10:ppn=1
#PBS -l walltime=50:00:00
#PBS -l mem=1200GB
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
filter_file="ext/hg19_filter_regions.txt"
steps="hic interaction postprocess"

mpirun -np 10 python ./snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps --parallel
