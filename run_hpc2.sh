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

indir=""
suffix=""
outdir=""
chrs=""
pos=""
chrlen="ext/mm10.chrom.sizes"
genome="mm10"
filter_file="ext/mm10_filter_regions.txt"
steps="hic interaction postprocess"

mpirun -np 10 python ./snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps --parallel

