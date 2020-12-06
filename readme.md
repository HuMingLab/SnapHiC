# SnapHiC: Single Nucleus Analysis Pipeline for Hi-C Data 
## Identifying chomatin loops from single cell Hi-C data
### 1. Installation 
1. Install python version >= 3.6 
2. Make sure MPI (e.g. [open-mpi](https://www.open-mpi.org/)) is installed on your system and is on your path (`which mpicc` should return a path). 
3. Use pip to install the required modules: 
```
pip install -r requirements.txt
```

### 2. Required Input Files:
1. A "tab-separated" or "tab-separated and gzipped" file of the reads mapped on the genome for each cell. (You can generate these files by running bwa-mem on your fastq files).
2. chrom.sizes file for the genome of interest. (can be downloaded from [here](https://hgdownload.soe.ucsc.edu/downloads.html)). Files for mm10 and hg19 are included in the `ext` directory. 
3. Optionally a binned bed file of filtered regions for the genome (aka blacklist regions). Filtered regions and low mappability regions for hg19 and mm10 in 10KB resolution are included in the `ext` directory. 

### 3. Run
We strongly recommend using an HPC environmet where you can request nodes/processors and allocate memory. However, runfiles for multi-threaded runs and single-processor systems are also provided.
1. Put all the mapped read files in one directory: one file for each cell, each file containing one line per read pair. For each read pair, the file should contain at least 2 columns for chromsomes and 2 columns for the mapped positions (bp).  
2. Either one of the following options will work but the HPC cluster is the recommended method:  
&Tab;**(I) For HPC clusters with a job scheduler such as PBS or SLURM**:  
&Tab;&Tab; Use the template provided for the PBS scheduler: *run_hpc1.sh*, followed by *run_hpc2.sh*. Set the variables in these files based on the point (3) below. (These files use the --parallel flag)  
&Tab;**(II) For a compute node with multiple cores but no scheduler**:  
&Tab;&Tab; Use the *run_threaded.sh* file as the template. Set the variables in this file based on the point (3) below. (This file uses --threaded flag, along with the number of threads you specify to use).  
&Tab;**(III) For a single core run non-parallel run**:  
&Tab;&Tab; This can be extremely slow. We do not recommended to use this mode for a large number of cells. Use *run_desktop.sh*.
3. Regardless of which running method you are using, the following variables have to be set in the corresponding run files:  
&Tab;`indir`="/path/where/the/mapped/data/are/stored"   
&Tab;`suffix`="contacts.txt" (Filename suffix for mapped files. This is used to distinguish input files if there are other files in the input directory)  
&Tab;`outdir`="/path/where/the/output/will/be/saved"  
&Tab;`chrs`="2 4" (Two integers indicating the column numbers of the chromosomes in the input files. Starting from 1).  
&Tab;`pos`="3 5" (Two integers indicating the column numbers of the read mapped positions in the input files. Starting from 1). 
&Tab;`chrlen`="ext/mm10.chrom.sizes" (chrom.sizes file. You can download it from the UCSC webpage if needed).  
&Tab;`genome`="mm10" (Name of the reference genome. Currently accepts mm10 and hg19). It is used to determine name of the chromosomes to process, as well as to create juicebox''s hic file.   
&Tab;`filter_file`="ext/mm10_filter_regions.txt" (this is optional. We recommend providing a binned bed file for areas of low mappability or filtered regions of the genome. We have included these files for mm10 and hg19 in the *ext* directory) 
&Tab;Additionally for the threaded run (single node with no scheduler) you will need to specify *num_proc* - number of threads to use).  
4. Execute the run file with modified variables (or submit it to the job scheduler). 

### 4. The Output: 
The final detected significant interactions are stored in the file: *<outdir>/postprocessed/combined.postprocessed.summits.bedpe*. This tab-separated file includes the following columns:  
- chr1, x1, x2, chr2, y1, y2: start and end position of interacting pair of regions. 
- outlier_count: number of cells with > 1.96 normalized rwr score  
- pvalue, tstat, fdr_dist: statistical measures computed for each binpair 
- case_avg, control_avg: average over cells signal of each binpair (case) and average signal of its negihboring region (control). 
- circle, donut, horizontal, vertical, lower_left: extra filters computed for each binpair. See the manuscript for more details. 
- i, j, min_dist, ro, rownum, delta, ref_neighbor, eta, rank, transformed_rank, transformed_eta: extra columns used in computation of the interaction clusters. 
- eta_cluster: cluster number (ID), IDs are chromosome specific. 
- cluster_size: number of interacting binpairs within a cluster. 
- neg_log10_fdr: measure of the strength of the cluster.  This number is the same for all binpairs within a cluster. We recommend using *fdr_dist* as a measure of strength of the interaction as opposed to this column. 
- summit: all binpairs in the summits file have a value of 1. Note that there might be multiple summits per cluster. 

In addition to the final output file metioned above, there are multiple other intermediate output files. Each of the steps in the process (bin, rwr, hic, interaction, postprocess) creates a separate subdirectory in the *<outdir>*. The first two steps generate output files per cell, and per cell-chromosome respectively. The remaining steps use the combined data and generate one file per chromosome. 

### 5. Testing  
You can use the dataset provided [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC_100.tar.gz) to test your installation of the program. This dataset contains a subset of 100 cells from the Ecker''s ODC data (hg19). The output for this set can be downloaded from [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/output/Ecker/ODC_100_output.tar). 
After downloading the input file, untar it so that you can see the 100 input data files. Run files that we have used to generate these results are also included with the snapHiC package in *run_test_odc_step1.sh* and *run_test_odc_step2.sh* scripts. Set the input/output directories and submit them for run.

### 6. Recommendations for parallel setting:  
As mentioned earlier, you would want to use as many processors as possible for the RWR step (first script), but you want to provide enough memory for the rest of the steps.  
The specific size of memory and processors per node will be different based on the genome being used and the number of cells in the dataset. Here are the setting we have used for our datasets.  
| genome | #cells | binsize | distance | run step | nodes | processors per node | memory per node | runtime |  
| --- | --- | --- | --- | --- | --- | --- | --- | --- |  
| mm10 | 100 | 10,000 | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 2.4 hrs |  
| mm10 | 100 | 10,000 | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 0.7 hrs |  
| mm10 | 200 | 10,000 | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 8.8 hrs |  
| mm10 | 200 | 10,000 | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 1.3 hrs |  
| mm10 | 300 | 10,000 | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 13.2 hrs |  
| mm10 | 300 | 10,000 | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 1.8 hrs |  
| mm10 | 400 | 10,000 | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 18.4 hrs |  
| mm10 | 400 | 10,000 | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 2.4 hrs |  
 
### 7. More Details on running:
In brief, you will need to run the program in the following two steps. These are the already included in the run-files, and by executing them as described above or by submitting to the job-scheduler you are essentially following these steps:
step1:
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column-numbers-for-chromosomes-in-data> -p <column-numbers-for-read-positions-in-data> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "bin rwr"
```
followed by step2:
```
python3 snap.py -i <input-directory> -o <output-directory> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "hic interaction postprocess"
```

The operation in this pipeline is divided into 5 steps: (1) binning, (2) RWR computation, (3) combining cells, (4) computation of local backgroud, (5) finalizing and postprocessing. You can specify which steps you want to run as a command line argument. To specify the steps use the argumet `--step` followed by any of the 'bin', 'rwr', 'hic', 'interaction', and 'postprocess'. (if you don''t enter any steps, it will attempt to run the entire pipeline). 
We strongly recommend breaking the operation down into two parts: the first part is the binning and RWR computation and the second part includes all the remaining steps. 
This breakdown is recommended for two reasons. First, when you have a large number of cells to process, you would want to run the operation in parallel as much as possible, increasing the number of processors per node; but the available memory is limited and your run might fail due to out of memory error (which the program is not capable of catching). By waiting until the RWR step is completed you can make sure all the cells and chromosomes are processed. **If you restart this rwr step, it will recognize the chromosomes/cells that have already been processed and won''t repeat the computation for them.** In case you run out of memory and you suspect some of the RWR computations might be corrupted, we have included a script (utils/validate_rwr.py) that can be used to validate that all RWRs are computed successfully. 
Second, for steps 3 and 4 you will need more memory per CPU, as such if your dataset contains more than 100 cells, we recommend using 1, 2, or 3 processors per node for a node with ~120GB memory. By breaking the operation down into these two steps, you can optimally adjust your memory and processor per node according to your HPC and the number of cells in your dataset.  

Additional arguments that can be used, can be seen by calling *python snap.py --help*. 

For additional questions, please submit an issue providing the details of your system and run. You can contact us via email <hum@ccf.org> (Ming Hu) or <a.abnousi@gmail.com> (Armen Abnousi).
