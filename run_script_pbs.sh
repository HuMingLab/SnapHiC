#!/bin/bash
#PBS -q home
#PBS -N test
#PBS -l nodes=2:ppn=16
#PBS -l walltime=30:00:00
#PBS -l pmem=5GB
#PBS -o log.${PBS_JOBID}.log
#PBS -e log.${PBS_JOBID}.log.err
#PBS -V
#PBS -M email@email.com
#PBS -m abe

cd /oasis/tscc/scratch/abnousa/snapHiC

num_processor=64
indir="/oasis/tscc/scratch/abnousa/snapHiC/inputs/F123_unphased"
suffix="rm_hotspot.sorted.txt"
outdir="/oasis/tscc/scratch/abnousa/snapHiC/output/F123_unphased"
chrs="3 7"
pos="4 8"
chrlen="/oasis/tscc/scratch/abnousa/snapHiC/inputs/chrlens/mm10.chrom.sizes"
genome="mouse"
#chrom="chr21"
dist=1000000
bin=10000
local_neighborhood_lower=3
local_neighborhood_upper=5
fdr_thresh=0.01
postproc_gap_large=5
postproc_gap_small=2
candid_lower_dist=100000
candid_upper_dist=900000
cluster_gap=3
max_memory=1

mpirun -np $num_processor python3 ./snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --dist $dist --binsize $bin --local-lower-limit $local_neighborhood_lower --local-upper-limit $local_neighborhood_upper --fdr-threshold $fdr_thresh --postproc-gap-large $postproc_gap_large --postproc-gap-small $postproc_gap_small --candidate-lower-distance $candid_lower_dist --candidate-upper-distance $candid_upper_dist --clustering-gap $cluster_gap --max-memory $max_memory --parallel
