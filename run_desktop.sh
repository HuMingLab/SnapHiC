#!/bin/bash

indir=""
suffix=""
outdir=""
chrs=""
pos=""
chrlen="ext/mm10.chrom.sizes"
genome="mm10"
fdr_thresh=0.1
filter_file="ext/mm10_filter_regions.txt"
steps="bin rwr hic interaction postprocess"

python ./snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --fdr-threshold $fdr_thresh --filter-file $filter_file --steps $steps
