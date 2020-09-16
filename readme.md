# snapHiC 
## Chromatin Interaction Caller For Single Cell HiC Data
### 1. Installation
Prior to run, please ensure following modules are installed in your python3 environment:
- numpy
- scipy
- statsmodels
- pandas
- skimage (scikit-image)
- networkx
- mpi4py (optional but recommended; required for run in HPC environment)

### 2. Run
#### TLDR;
1. Put all the mapped files in one directory. For each read pair your file should have at least 2 columns for chromsomes and 2 columns for the mapped position (bp).  
2. Either one of the following options will work but the HPC cluster is the recommended method:  
&Tab;**(I) For HPC clusters with a job scheduler such as PBS or SLURM**:  
&Tab;&Tab; Use the template provided for the PBS scheduler: *run_hpc1.sh*, followed by *run_hpc2.sh*. Set the variables in these files based on the point (3) below. (These files use the --parallel flag)  
&Tab;**(II) For a compute node with multiple cores but no scheduler**:  
&Tab;&Tab; Use the *run_threaded.sh* file as template. Set the variables in this file based on the point (3) below. (This file uses --threaded flag, along with the number of threads you specify to use).  
&Tab;**(III) For a single core run non-parallel run**:  
&Tab;&Tab; This can be extremely slow, it is not recommended to use this mode for large number of cells. Use *run_desktop.sh*.
3. Regardless of which running method you are using, the following variables have to be set in the corresponding run files:  
&Tab;`indir`="/path/where/the/mapped/data/are/stored"   
&Tab;`suffix`="contacts.txt" (filename suffix for mapped files. This is used to distinguish input files if there are other files in the input directory)  
&Tab;`outdir`="/path/where/the/output/will/be/saved"  
&Tab;`chrs`="2 4" (two integers indicating column numbers of the chromosomes in the input files. starting from 1).  
&Tab;`pos`="3 5" (two integers indicating column numbers of the read mapped locations in the input files. starting from 1).  
&Tab;`chrlen`="ext/mm10.chrom.sizes" (chrom.sizes file. You can download from ucsc webpage if needed).  
&Tab;`genome`="mm10" (name of the chromosome. Currently accepts mmxx and hgxx. Is used to determine name of chromosomes to process)  
&Tab;`fdr_thresh`=0.1 (depeding on your dataset we recommend using 0.1 or 0.01 threshold for FDR)  
&Tab;`filter_file`="ext/mm10_filter_regions.txt" (this is optional. We recommend providing a binned bed file for areas of low mappability or filtered regions of the genome. We have included these files for mm10 and hg19 in the *ext* directory)  
&Tab;Additionally for the threaded run (single node with no scheduler) you will need to specify *num_proc* - number of threads to use).  

#### Details on running:  
In brief, you will need to run the program in the following two steps:   
step1:  
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column-numbers-for-chromosomes-in-data> -p <column-numbers-for-read-positions-in-data> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "bin rwr"
```
followed by step2:  
```
python3 snap.py -i <input-directory> -o <output-directory> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "hic interaction postprocess"
```

The operation in this pipeline is divided into 5 steps: (1) binning, (2) RWR computation, (3) combining cells, (4) computation of local backgroud, (5) finalizing and postprocessing. You can specify which steps you want to run as a command line argument (if you don''t enter any steps, it will attempt to run the entire pipeline).  
However, we strongly recommend breaking the operation down into two parts, in the first part is the binning and RWR computation and in the second part the remaining steps. This is recommended for two reasons. First, when you have a large number of cells to process, you would want to run the operation in parallel as much as possible, increasing the number of processors per node; but the available memory is limited and your run might fail due to out of memory error (which the program is not capable of catching). By waiting until the RWR step is completed you can make sure all the cells and chromosomes are processed. If you restart this step, it will recognize the chromosome/cells that have already been processed and won''t repeat the computation for them. In case you run out of memory and you suspect some of the RWR computations might be corrupted, we have included a script (utils/validate_rwr.py) that can be used to validate that all RWRs are computed successfully. Second, for steps 3 and 4 you will need more memory per CPU, as such if your dataset contains more than 100 cells, we recommend using 1 or 2 processor per node for a node with ~120GB memory. By breaking the operation down into these two steps, you can optimally adjust your memory and processor per node according to your HPC and the number of cells in your dataset.  

Additional arguments that can be used, can be seen by calling *python snap.py --help*.
Use the following command to run the software:
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column numbers for chromosomes in data> -p <column numbers for read positions in data> -l <chromosome lengths file> -g <genome: human or mouse>
```
Please note that in order to run the code in parallel mode (e.g. in cluster computing systems), you will have to use the `--parallel` argument when calling the software.   
A complete listing of all the permitted parameters to the software can be seen by using `--help` argument.  

### 3. Testing  
You can use the dataset provided [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC) to test your installation of the program. This dataset contains a subset of 100 cells from the Ecker''s ODC data (hg19). The output for this set is [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC).   

### 4. Recommendations for parallel setting:  
As mentioned earlier, you would want to use as many processors as possible for the RWR step (first script), but you want to provide enough memory for the rest of the steps.  
The specific size of memory and processors per node will be different based on the genome being used and the number of cells in the dataset. Here are the setting we have used for our datasets.  
| genome | #cells | binsize | run step | nodes | processors per node | memory per node | runtime |  
| --- | --- | --- | --- | --- | --- | --- | --- |  
| hg19 | 100 | 10,000 | 1 (bin rwr) | 15 | 2 | 120GB | 3 hrs |  
| hg19 | 100 | 10,000 | 2 (hic inter. postproc.) | 10 | 1 | 120GB | 4 hrs |  
| mm10 | 100 | 10,000 | 1 (bin rwr) | 10 | 3 | 120GB | 2 hrs |  
| mm10 | 100 | 10,000 | 2 (hic inter. postproc.) | 10 | 2 | 120GB | 2 hrs |

