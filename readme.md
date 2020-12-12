# SnapHiC: Single Nucleus Analysis Pipeline for Hi-C Data 
#### (Latest updates: December 7th, 2020, version 0.1.0)
## Identifying chromatin loops from single cell Hi-C data
### 1. Installation
You can download the singularity image/recipe of SnapHiC [here](http://renlab.sdsc.edu/abnousa/snapHiC/singularity_releases). If you don't have access to singularity or prefer your own installation, please follow the following steps:    
1. Install python version >= 3.6 
2. Make sure MPI (e.g. [open-mpi](https://www.open-mpi.org/)) is installed on your system and is on your path (`which mpicc` can return a path). 
3. Use pip to install the required modules: 
```
pip install -r requirements.txt
```
4. Install Java version >= 8 and include it in your path. Java is required to generate .hic files for visualization in Juicebox.

### 2. Required input files:
1. "Tab-separated" or "tab-separated and gzipped" files containing the mapped read pairs (contacts) for each single cell. These contact files can be generated from raw fastq files following the methods decribed in previous publications ([PMID: 31384045](https://pubmed.ncbi.nlm.nih.gov/31384045/), [PMID: 31501549](https://pubmed.ncbi.nlm.nih.gov/31501549/), and [PMID: 28682332](https://pubmed.ncbi.nlm.nih.gov/28682332/)), or other single cell Hi-C data preprocessing pipelines such as Dip-C (https://github.com/tanlongzhi/dip-c). In each file, one line represents one contact with 2 columns for chromosome name and 2 columns for the mapped positions (bp). 
2. chrom.sizes file for the genome of interest, which can be downloaded from [here](https://hgdownload.soe.ucsc.edu/downloads.html). Files for mm10, hg19 and hg38 are included in the `ext` directory. 
3. A binned bed file of the genomic regions that are excluded from loop calling. In our study, we defined mappability for each 10KB bin based on the restriction enzyme MboI, as described in our previous study [PMID: 23023982](https://pubmed.ncbi.nlm.nih.gov/23023982/). We removed all 10KB bins with mappability <=0.8, and all 10KB bins overlapped with the ENCODE blacklist regions. Filtered regions for mm10, hg19 and hg38 at 10KB resolution with the restriction enzyme MboI, i.e., bins with mappability <=0.8 and bins overlapped with the ENCODE blacklist regions, are included in the `ext` directory. Local genomic features (including mappability scores) for different reference genomes, different bin resolutions, and different restriction enzymes can be downloaded [here](http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/).

### 3. Running SnapHiC
We strongly recommend using an HPC environment where you can request multiple nodes/processors and allocate memory. Alternatively, you can still run snapHiC using multithreaded and single-processor environments, provided enough memory and runtime.  
1. Put all the contact files, one for each cell, into the same directory.  
2. Edit the *run_step1.sh* and *run_step2.sh* files. If you use an HPC cluster with a job scheduler such as PBS or SLURM (we strongly recommend), specify the required nodes, processors, memory, and load the required modules (python3.6+, MPI, java8+). If you use a regular compute node without job scheduler, or a desktop computer, you can skip this step (you still need big memory and the computation can be slow).    
3. Set the following variables in the *run_step1.sh* and *run_step2.sh* files:  
&Tab;`snapHiC_dir`="/path/to/directory/where/snapHiC/is/located/" (path to the directory contains the *snap.py* file in SnapHiC).  
&Tab;`parallelism`="parallel" (it can take one of the three options: **parallel**, **threaded**, or **single-proc**. Use **parallel** if you use an HPC with job scheduler, **threaded** if you use multiple processors without job scheduler, and **singl-proc** otherwise).    
&Tab;`number_of_processors`=15 (if you use **threaded** or **parallel**, please specify the number of processors).  
&Tab;`indir`="/path/where/the/contact/files/are/stored" (contact files should be tab separated. They can be gzipped).  
&Tab;`suffix`="contacts.txt.gz" (File name suffix for the contact files, which is used to distinguish input files if there are other files in the same input directory).  
&Tab;`outdir`="/path/where/the/output/will/be/saved"  
&Tab;`chrs`="2 4" (Two integers indicating the column numbers of the chromosomes in the contact files. Starting from 1).  
&Tab;`pos`="3 5" (Two integers indicating the column numbers of the read mapped positions in the contact files. Starting from 1).  
&Tab;`chrlen`="ext/mm10.chrom.sizes" (chrom.sizes file. You can download it from [here](https://hgdownload.soe.ucsc.edu/downloads.html)).  
&Tab;`genome`="mm10" (Name of the reference genome. SnapHiC currently accepts mm10, hg19 and hg38). It is used to determine the number of autosomal chromosomes, and to generate .hic file for visualization in Juicebox.   
&Tab;`filter_file`="ext/mm10_filter_regions.txt" (A binned bed file of the genomic regions that are excluded from loop calling, such as the low mappability regions and the ENCODE blacklist regions. We provide these files for mm10, hg19 and hg38 at 10KB resolution with the restriction enzyme MboI in the *ext* directory).   
&Tab;If you use **threaded** (single node without job scheduler), you need to specify *num_proc*, the number of threads to use.  
4. Execute the run file with the modified variables, or submit it to the job scheduler. 

### 4. The output file: 
The SnapHiC-identified chromatin loops are stored in the file: *<outdir>/postprocessed/*.postprocessed.summits.bedpe*. (* will be replaced by the library name provided as input via the **--prefix** argument. If no such argument is provided, it will be replaced with *combined*) This tab-separated file includes the following columns:  
- chr1, x1, x2, chr2, y1, y2: start and end position of chromatin loops. 
- outlier_count: the number of cells with >1.96 normalized contact probability (with respect to global background) at the loop summit.  
- pvalue, tstat, fdr_dist: statistical measures (P-values from the paired t-test, t-statistics from the paired t-test and false discovery rate for all bin pairs at the same 1D genomic distance, with respect to local background) computed for each chromatin loop. We recommend using fdr_dist as the measure of the chromatin loop strength.
- case_avg, control_avg: across all cells, the average normalized contact probability of each chromatin loop (case) and its local neighboring region (control). 
- circle, donut, horizontal, vertical, lower_left: the average number of cells with >1.96 normalized contact probability at the five local background regions. SnapHiC applies extra folder change filters with respect to five local background regions. See details in the manuscript. 
 
Additionaly, a second file containing all the candidates (less sigletons), will be stored in: *<outdir>/postprocessed/*.postprocessed.all_candidates.bedpe*. In addition to the columns above, this file will contain the following extra columns:  
- eta_cluster: cluster number (ID), which is chromosome specific. 
- cluster_size: the number of loop candidates within a cluster. 
- neg_log10_fdr: measure of the strength of the cluster. This value is the same for all loop candidates within a cluster. We recommend using *fdr_dist* as the measure of the chromatin loop strength. 

SnapHiC also generates .hic and .cool files stored at *<outdir>/hic/allChr.[hic/cool]*. SnapHiC first computes the % of outlier cells (i.e., the proportion of cells with normalized contact probability >1.96), and then takes the integer ceiling of 100*(% of outlier cells) to create a count matrix. SnapHiC then uses the Juicer software and the cooler software to convert the count matrix into .hic file and .cool file, respectively. User can skip this step by including **--no-hic** or **--no-cool** argument in the execution. 

In addition to the output files mentioned above, SnapHiC generates multiple intermediate output files. Each of the five steps in the process (bin, rwr, hic, interaction and postprocess) creates a separate sub-directory in the `outdir`. The first two steps (bin and rwr) generate output files for each chromosome in each cell. The remaining three steps (hic, interaction and postprocess) combine all cells and generate one file for each chromosome. 

### 5. Testing SnapHiC 
To test SnapHiC, please download this sample dataset [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC_100.tar.gz) (1.2GB), which contains the contact files of 100 randomly selected oligodendrocytes from Lee et al study ([PMID: 31501549](https://pubmed.ncbi.nlm.nih.gov/31501549/), Ref: hg19). The final list of 6,249 loops can be downloaded from [here](http://renlab.sdsc.edu/hum/ODC_100_summits.bedpe) (2.4MB). The complete output files for this sample dataset can be downloaded from [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/output/Ecker/ODC_100_output.tar) (451GB). 

After downloading this sample dataset, untar it to generate 100 contact files, each containing the contacts for one cell. Run-files for this sample dataset are in the *test_run_scripts* directory. You can set the input and output directories, load the required modules, and submit the run-files.

### 6. Recommendations for parallel setting:  
You can use as many processors as possible for the RWR step, as long as each processor has sufficient memory. For 10Kb bin resolution, we recommend 30-40GB memory per processor. The remaining steps are performed for each chromosome separately, therefore, we recommend using as many processors as the number of chromosomes in the genome. Providing sufficient memory for each processor for the remaining steps is important. For 2Mb maximal genomic distance, 10Kb bin resolution, and for human and mouse genome, we recommend providing 50GB memory per processor. Here are the settings we have used in our study.  
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
 
### 7. More details on running SnapHiC:
The operation in SnapHiC consists of five steps: (1) binning, (2) random walk with restart (RWR) computation, (3) combining cells, (4) computation of local background, (5) finalizing and postprocessing. You can specify which steps you want to run as a command line argument. To specify the steps, you can use the argument `--step` followed by any of the five options: 'bin', 'rwr', 'hic', 'interaction', and 'postprocess'. If you don't specify any steps, SnapHiC will run all five steps. 

We strongly recommend running SnapHiC in the following two parts, which are described in the run-files. 

Part 1:
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column-numbers-for-chromosomes-in-data> -p <column-numbers-for-read-positions-in-data> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "bin rwr"
```

followed by Part 2:
```
python3 snap.py -i <input-directory> -o <output-directory> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "hic interaction postprocess"
```

First of all, when your dataset consists of a large number of cells, you can run the computation in parallel as much as possible by using a large number of processors per node. However, since the available memory is limited, your computation might fail due to the out of memory error (which SnapHiC currently is not able to detect). By waiting until the RWR step is finished, you can make sure that all chromosomes in all cells are processed. **If you restart the RWR step, SnapHiC can recognize the chromosomes/cells that have already been processed and will not repeat the RWR computation for them.** 

In case you run out of memory and you suspect that some of the RWR computations might be corrupted, you can use the script (utils/validate_rwr.py) to assess whether RWR computation finished successfully in all cells. The input of the **validate_rwr.py** script is the output directory of the SnapHiC run. Please use command: `python utils/validate_rwr.py /path/to/output/of/your/run`. This script will generate a file named **missing.txt** in the <output>/rwr directory, which contains the name of RWR chromosomes/cells in which the RWR computation is corrupted. You can delete these files and restart the RWR step in SnapHiC.

In addition, for steps 3 and 4, you need more memory per processors. If your dataset contains more than 100 cells, we recommend using 1, 2, or 3 processors per node for a node with ~120GB memory. By breaking down the computation into these two parts, you can optimally adjust your memory and processor per node according to your HPC resource and the number of cells in your dataset.  

Additional arguments that can be used, which can be found by *python snap.py --help*. 

For any questions, comments and suggestions regarding SnapHiC, please submit an issue to Github with the details of your system and run. You can also send email to Armen Abnousi (<a.abnousi@gmail.com>) or Ming Hu (<hum@ccf.org>).
