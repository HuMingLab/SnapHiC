# snapHiC 
## Chromatin Interaction Caller For Single Cell HiC Data

Prior to run, please ensure following modules are installed in your python3 environment:
- numpy
- scipy
- statsmodels
- pandas
- skimage (scikit-image)
- networkx
- mpi4py (optional; required for run in HPC environment)

Use the following command to run the software:
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column numbers for chromosomes in data> -p <column numbers for read positions in data> -l <chromosome lengths file> -g <genome: human or mouse>
```
Please note that in order to run the code in parallel mode (e.g. in cluster computing systems), you will have to use the `--parallel` argument when calling the software.  
A complete listing of all the permitted parameters to the software can be seen by using `--help` argument.
```
python3 snap.py --help
usage: snap.py [-h] -i INDIR -s SUFFIX -o OUTDIR -c CHR_COLUMNS CHR_COLUMNS -p
               POS_COLUMNS POS_COLUMNS -l CHR_LENS [-g GENOME] [--chrom CHROM]
               [--dist DIST] [--binsize BINSIZE] [--low-cutoff LOW_CUTOFF]
               [--alpha ALPHA] [--parallel] [--threaded] [-n NUM_PROC]
               [--outlier OUTLIER] [--case-control-diff CASE_CONTROL_DIFF]
               [--local-lower-limit LOCAL_LOWER_LIMIT]
               [--local-upper-limit LOCAL_UPPER_LIMIT]
               [--fdr-threshold FDR_THRESHOLD]
               [--postproc-gap-large POSTPROC_GAP_LARGE]
               [--postproc-gap-small POSTPROC_GAP_SMALL]
               [--candidate-lower-distance CANDIDATE_LOWER_DISTANCE]
               [--candidate-upper-distance CANDIDATE_UPPER_DISTANCE]
               [--circle-threshold-multiplier CIRCLE_THRESHOLD_MULTIPLIER]
               [--donut-threshold-multiplier DONUT_THRESHOLD_MULTIPLIER]
               [--lower-left-threshold-multiplier LOWER_LEFT_THRESHOLD_MULTIPLIER]
               [--horizontal-threshold-multiplier HORIZONTAL_THRESHOLD_MULTIPLIER]
               [--vertical-threshold-multiplier VERTICAL_THRESHOLD_MULTIPLIER]
               [--outlier-threshold-multiplier OUTLIER_THRESHOLD_MULTIPLIER]
               [--summit-gap SUMMIT_GAP] [--filter-file FILTER_FILE]
               [--clustering-gap CLUSTERING_GAP] [--max-memory MAX_MEMORY]
               [--steps [STEPS [STEPS ...]]]

optional arguments:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        input directory
  -s SUFFIX, --suffix SUFFIX
                        suffix of the input files
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -c CHR_COLUMNS CHR_COLUMNS, --chr-columns CHR_COLUMNS CHR_COLUMNS
                        two integer column numbers for chromosomes
  -p POS_COLUMNS POS_COLUMNS, --pos-columns POS_COLUMNS POS_COLUMNS
                        two integer column numbers for read positions
  -l CHR_LENS, --chr-lens CHR_LENS
                        path to the chromosome lengths file
  -g GENOME, --genome GENOME
                        genome name; mouse/human
  --chrom CHROM         chromosome to process
  --dist DIST           distance from diagonal to consider
  --binsize BINSIZE     bin size used for binning the reads
  --low-cutoff LOW_CUTOFF
                        cut-off for removing short-range reads
  --alpha ALPHA         restart probability of random walk
  --parallel            if set, will attempt to run in parallel mode
  --threaded            if set, will attempt to use multiprocessing on single
                        machine
  -n NUM_PROC, --num-proc NUM_PROC
                        number of processes used in threaded mode
  --outlier OUTLIER     percentage threshold for finding outliers.
  --case-control-diff CASE_CONTROL_DIFF
                        difference threshold for case mean to control mean for
                        finding candidates
  --local-lower-limit LOCAL_LOWER_LIMIT
                        number of bins around center (in each direction) to
                        exlude from neighborhood
  --local-upper-limit LOCAL_UPPER_LIMIT
                        number of bins around center (in each direction)
                        forming the neighborhood
  --fdr-threshold FDR_THRESHOLD
                        FDR threshold used for candidate peak detection
  --postproc-gap-large POSTPROC_GAP_LARGE
                        number of bins around peaks to consider in
                        postprocessing
  --postproc-gap-small POSTPROC_GAP_SMALL
                        number of bins around peaks to disregard in
                        postprocessing
  --candidate-lower-distance CANDIDATE_LOWER_DISTANCE
                        lower threshold for distance between candidate peak
                        binpairs
  --candidate-upper-distance CANDIDATE_UPPER_DISTANCE
                        upper threshold for distance between candidate peak
                        binpairs
  --circle-threshold-multiplier CIRCLE_THRESHOLD_MULTIPLIER
                        multiplier for circle filter threshold used in finding
                        candidates
  --donut-threshold-multiplier DONUT_THRESHOLD_MULTIPLIER
                        multiplier for donut filter threshold used in finding
                        candidates
  --lower-left-threshold-multiplier LOWER_LEFT_THRESHOLD_MULTIPLIER
                        multiplier for lowerleft filter threshold used in
                        finding candidates
  --horizontal-threshold-multiplier HORIZONTAL_THRESHOLD_MULTIPLIER
                        multiplier for horizontal filter threshold used in
                        finding candidates
  --vertical-threshold-multiplier VERTICAL_THRESHOLD_MULTIPLIER
                        multiplier for vertical filter threshold used in
                        finding candidates
  --outlier-threshold-multiplier OUTLIER_THRESHOLD_MULTIPLIER
                        multiplier for number of outlier cells for threshold
                        used in finding candidates
  --summit-gap SUMMIT_GAP
                        disallowed distance between summit points of a cluster
  --filter-file FILTER_FILE
                        bed file of regions to be filtered. Regions should be
                        binned
  --clustering-gap CLUSTERING_GAP
                        number of allowed gaps between peaks in same cluster
  --max-memory MAX_MEMORY
                        memory available in GB, that will be used in
                        constructing dense matrices
  --steps [STEPS [STEPS ...]]
                        steps to run. Combination of bin,rwr,hic,interaction,
                        and postprocess are accepted (no comma, space
                        separated). Default is all steps.
```
