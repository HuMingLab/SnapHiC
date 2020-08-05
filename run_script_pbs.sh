#!/bin/bash
#PBS -q home
#PBS -N test
#PBS -l nodes=2:ppn=8
#PBS -l walltime=30:00:00
#PBS -l mem=240GB
#PBS -o log.${PBS_JOBID}.log
#PBS -e log.${PBS_JOBID}.log.err
#PBS -V
#PBS -M email@email.com
#PBS -m abe

cd /oasis/tscc/scratch/abnousa/snapHiC

num_processor=16
indir="/oasis/tscc/scratch/abnousa/snapHiC/inputs/F123_unphased"
suffix="rm_hotspot.sorted.txt"
outdir="/oasis/tscc/scratch/abnousa/snapHiC/output/F123_unphased"
chrs="3 7"
pos="4 8"
chrlen="/oasis/tscc/scratch/abnousa/snapHiC/inputs/chrlens/mm10.chrom.sizes"
genome="mouse"
genomeID="mm10"
chrom="None" #None or no argument for all chromosomes
dist=1000000
bin=10000
local_neighborhood_lower=2
local_neighborhood_upper=5
outlier_percent=0.97500211 # ~=1.96
case_control_diff=0.4
fdr_thresh=0.01
postproc_gap_large=5
postproc_gap_small=2
candid_lower_dist=100000
candid_upper_dist=900000
circle_threshold_multiplier=1.33
donut_threshold_multiplier=1.33
lower_left_threshold_multiplier=1.33
horizontal_threshold_multiplier=1.2
vertical_threshold_multiplier=1.2
outlier_threshold_multiplier=0.1
summit_gap=10000
filter_file="ext/mm10_filter_regions.txt"
cluster_gap=3
max_memory=1
steps="bin rwr hic interaction postprocess"

mpirun -np $num_processor python3 ./snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --dist $dist --binsize $bin --local-lower-limit $local_neighborhood_lower --local-upper-limit $local_neighborhood_upper --outlier $outlier_percent --case-control-diff $case_control_diff --fdr-threshold $fdr_thresh --postproc-gap-large $postproc_gap_large --postproc-gap-small $postproc_gap_small --candidate-lower-distance $candid_lower_dist --candidate-upper-distance $candid_upper_dist --circle-threshold-multiplier $circle_threshold_multiplier --donut-threshold-multiplier $donut_threshold_multiplier --lower-left-threshold-multiplier $lower_left_threshold_multiplier --horizontal-threshold-multiplier $horizontal_threshold_multiplier --vertical-threshold-multiplier $vertical_threshold_multiplier --outlier-threshold-multiplier $outlier_threshold_multiplier --summit-gap $summit_gap --filter-file $filter_file --clustering-gap $cluster_gap --max-memory $max_memory --chrom $chrom --steps $steps --parallel 

############################################
infile=$outdir"/hic/allChr.hic.input"
outfile=$outdir"/hic/allChr.hic"
java -jar ./utils/juicer_tools_1.22.01.jar pre $infile $outfile $genomeID
