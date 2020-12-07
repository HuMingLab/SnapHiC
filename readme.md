# SnapHiC: Single Nucleus Analysis Pipeline for Hi-C Data (Latest updates: December 7th, 2020)
## Identifying chomatin loops from single cell Hi-C data
### 1. Installation
You can download the singularity image of SnapHiC [here](http://renlab.sdsc.edu/abnousa/snapHiC/singularity_releases) to avoid installation. If you don't have access to singularity or prefer your own installation, please follow the following steps:    
1. Install python version >= 3.6 
2. Make sure MPI (e.g. [open-mpi](https://www.open-mpi.org/)) is installed on your system and is on your path (`which mpicc` can return a path). 
3. Use pip to install the required modules: 
```
pip install -r requirements.txt
```

### 2. Required input files:
1. A "tab-separated" or "tab-separated and gzipped" file of the mapped reads for each cell. For example, you can generate these files by running bwa-mem on your fastq files.
2. chrom.sizes file for the genome of interest, which can be downloaded from [here](https://hgdownload.soe.ucsc.edu/downloads.html). Files for mm10, hg19 and hg38 are included in the `ext` directory. 
3. A binned bed file of the filtered regions for the genome (the blacklist regions). Filtered regions and low mappability regions for mm10, hg19 and hg38 at 10KB resolution with the restriction enzyme MboI are included in the `ext` directory. 

### 3. Run
We strongly recommend using an HPC environment where you can request nodes/processors and allocate memory. As alternatives, we also provide run-files for multi-threaded run and single-processor system.
1. Put all the mapped read files, one for each cell, into the same directory. In each file, one line represents one mapped read pair. For each mapped read pair, the file should contain 2 columns for chromsomes and 2 columns for the mapped positions (bp).  
2. Edit the *run_step1.sh* and *run_step2.sh* files. If you use an HPC cluster with a job scheduler such as PBS or SLURM (we strongly recommend), modify the first few lines to set the required nodes, processors, memory, and load the required modules (python3.6+, MPI, and the packages installed using pip as described in the Installation section above). If you use a regular compute node with no job scheduler, or a desktop computer, you can skip this step (you still need big memory and the computing is slow).    
3. Set the following variables in the *run_step1.sh* and *run_step2.sh* files:  
&Tab;`snapHiC_dir`="/path/to/directory/where/snapHiC/is/located/" (path to the directory contains snap.py file of the SnapHiC pipeline).  
&Tab;`parallelism`="parallel" (it can be one of the **parallel**, **threaded**, or **single-proc**. Use **parallel** if you use an HPC with job scheduler, **threaded** if you use multiple processors without job scheduler, and **singl-proc** otherwise).    
&Tab;`number_of_processors`=15 (if you use **threaded** or **parallel**, please specify the number of processors).  
&Tab;`indir`="/path/where/the/mapped/data/are/stored" (files should be tab separated. They can be gzipped).  
&Tab;`suffix`="contacts.txt.gz" (Filename suffix for the mapped read files, which is used to distinguish input files if there are other files in the same input directory).  
&Tab;`outdir`="/path/where/the/output/will/be/saved"  
&Tab;`chrs`="2 4" (Two integers indicating the column numbers of the chromosomes in the mapped read files. Starting from 1).  
&Tab;`pos`="3 5" (Two integers indicating the column numbers of the read mapped positions in the mapped read files. Starting from 1).  
&Tab;`chrlen`="ext/mm10.chrom.sizes" (chrom.sizes file. You can download it from [here](https://hgdownload.soe.ucsc.edu/downloads.html)).  
&Tab;`genome`="mm10" (Name of the reference genome. SnapHiC currently accepts mm10, hg19 and hg38). It is used to determine the number of autosomal chromosomes, and to generate .hic file for visualization in Juicebox.   
&Tab;`filter_file`="ext/mm10_filter_regions.txt" (This is optional. We recommend providing a binned bed file for regions with low mappability or filtered regions of the genome, such as the ENCODE blacklist regions and the MHC locus. We provided these files for mm10, hg19 and hg38 at 10Kb resolution with the restriction enzyme MboI in the *ext* directory).   
&Tab;If you use **threaded** (single node without job scheduler), you need to specify *num_proc* - the number of threads to use).  
4. Execute the run file with the modified variables, or submit it to the job scheduler. 

### 4. The output file: 
The SnapHiC-identified chromatin loops are stored in the file: *<outdir>/postprocessed/combined.postprocessed.summits.bedpe*. This tab-separated file includes the following columns:  
- chr1, x1, x2, chr2, y1, y2: start and end position of chromatin loops. 
- outlier_count: number of cells with >1.96 normalized contact probability.  
- pvalue, tstat, fdr_dist: statistical measures (P-values from the paired t-test, t-statistics from the paired t-test and false discovery rate at the same 1D genomic distance) computed for each chromatin loop. 
- case_avg, control_avg: across all cells, the average normalized contact probability of each chromatin loop (case) and its local neighboring region (control). 
- circle, donut, horizontal, vertical, lower_left: extra filters computed for each chromatin loop. See details in the manuscript. 
- i, j, min_dist, ro, rownum, delta, ref_neighbor, eta, rank, transformed_rank, transformed_eta: extra columns used to group nearby loop candidates into loop clusters. See details in the manuscript.
- eta_cluster: cluster number (ID), which is chromosome specific. 
- cluster_size: the number of loop candidates within a cluster. 
- neg_log10_fdr: measure of the strength of the cluster. This number is the same for all loop candidates within a cluster. We recommend using *fdr_dist* as the measure of the chromatin loop strength. 
- summit: all chromatin loops in the summits file have a value of 1. Note that there might be multiple summits for the same cluster. 

In addition to the output file metioned above, SnapHiC generates multiple intermediate output files. Each of the five steps in the process (bin, rwr, hic, interaction and postprocess) creates a separate sub-directory in the *<outdir>*. The first two steps (bin and rwr) generate output files for each chromosome in each cell. The remaining three steps (hic, interaction and postprocess) combine all cells and generate one file for each chromosome. 

### 5. Testing SnapHiC 
You can use this sample dataset [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC_100.tar.gz) to test your installation of the SnapHiC pipeline. This sample dataset contains 100 oligodendrocytes from Lee et al study (PMID: 31501549, Ref: hg19). The output for this sample dataset can be downloaded from [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/output/Ecker/ODC_100_output.tar). 
After downloading the input file, untar it so that you can find 100 input files, each representing the mapped read pairs for one cell. Run-files to generate these results are included in the SnapHiC package at the *test_run_scripts* directory. You can set the input/output directories, load the required modules, and submit the run-files.

### 6. Recommendations for parallel setting:  
You can use as many processors as possible for the RWR step (first script), as long as you provide sufficient memory for the remaining steps.  
The size of memory and the number of processors per node are different for different reference genomes and different numbers of cells in the dataset. Here are the setting we have used in our study.  
| genome | #cells | binsize | distance | run step | nodes | processors per node | memory per node | runtime |  
| --- | --- | --- | --- | --- | --- | --- | --- | --- |  
| mm10 | 100 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 2.4 hrs |  
| mm10 | 100 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 0.7 hrs |  
| mm10 | 200 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 8.8 hrs |  
| mm10 | 200 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 1.3 hrs |  
| mm10 | 300 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 13.2 hrs |  
| mm10 | 300 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 1.8 hrs |  
| mm10 | 400 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | 18.4 hrs |  
| mm10 | 400 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | 2.4 hrs |  
 
### 7. More details on running:
You need to run SnapHiC in the following two parts, which are described in the run-files. 
Part 1:
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column-numbers-for-chromosomes-in-data> -p <column-numbers-for-read-positions-in-data> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "bin rwr"
```
followed by Part 2:
```
python3 snap.py -i <input-directory> -o <output-directory> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "hic interaction postprocess"
```

The operation in SnapHiC consists of 5 steps: (1) binning, (2) random walk with restart (RWR) computation, (3) combining cells, (4) computation of local backgroud, (5) finalizing and postprocessing. You can specify which steps you want to run as a command line argument. To specify the steps, you can use the argument `--step` followed by any of the 'bin', 'rwr', 'hic', 'interaction', and 'postprocess'. If you don't specify any steps, SnapHiC will run all 5 steps. 

We strongly recommend breaking down the computation into two parts: the first part includes binning and RWR computation and the second part includes the remaining three steps. 

First of all, when your dataset consists of a large number of cells, you can run the computation in parallel as much as possible by using a large number of processors per node. However, since the available memory is limited, your computation might fail due to the out of memory error (which SnapHiC currently is not able to detect). By waiting until the RWR step is finished, you can make sure that all chromosomes in all cells are processed. **If you restart the rwr step, SnapHiC can recognize the chromosomes/cells that have already been processed and does not repeat the rwr computation for them.** In case you run out of memory and you suspect that some of the RWR computations might be corrupted, you can use the script (utils/validate_rwr.py) to assess whether RWR are computed successfully in all cells. 

In addition, for steps 3 and 4, you need more memory per processors. If your dataset contains more than 100 cells, we recommend using 1, 2, or 3 processors per node for a node with ~120GB memory. By breaking down the computation into these two parts, you can optimally adjust your memory and processor per node according to your HPC resource and the number of cells in your dataset.  

Additional arguments that can be used, which can be found by *python snap.py --help*. 

For additional questions, comments and suggestions, please submit an issue with the details of your system and run. You can contact us via email <hum@ccf.org> (Ming Hu) or <a.abnousi@gmail.com> (Armen Abnousi).
