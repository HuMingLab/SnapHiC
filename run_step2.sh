#!/bin/bash
#PBS -q hotel
#PBS -N snHiC_s2
#PBS -l nodes=5:ppn=3
#PBS -l walltime=50:00:00
#PBS -l mem=100GB
#PBS -j oe
#PBS -o log.${PBS_JOBNAME}.${PBS_JOBID}.log

### SET THE DIRECTIVES ABOVE IF YOU ARE WORKING IN HPC WITH A JOB SCHEDULER
### AND LOAD THE REQUIRED MODULES IF NEEDED. REQUIRES PYTHON3.6+ WITH 
### MODULES TO BE INSTALLED USING "pip install -r requirements.txt"
#cd ${PBS_O_WORKDIR}
#module load python/3.6.1

## Required for cooler:
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

############################################################################
###                            User Variables                            ###
############################################################################
snapHiC_dir="/projects/ps-renlab/abnousa/snapHiC"	#where the snapHiC is located on your system
parallelism="parallel" 					#options are "parallel" "threaded" "singleproc"
number_of_processors=15   				#required only if parallelism is set to parallel or threaded
indir="" 						#directory containing input files (e.g. *.pairs files)
suffix="" 						#all input files should have the same suffix. it can be an empty string "", or ".txt"
outdir="" 						#directory where output files will be stored
chrs="" 						#2 integer numbers, the column number of chromosomes in the input files. (e.g. "3 5") starting from 1
pos="" 							#2 integer numbers, the column number of mapped position of read-pairs. (e.g.  "4 6") starting from 1
chrlen="/projects/ps-renlab/abnousa/snapHiC/ext/mm10.chrom.sizes" 		#path to the chrom.sizes file"
genome="mm10"  						#genomeID that will be used for genereation of ".hic" file 
filter_file="/projects/ps-renlab/abnousa/snapHiC/ext/mm10_filter_regions.txt" 	#regions to be filtered, for example due to low mappability
steps="hic interaction postprocess" 					#steps to run the pipeline. Recommended (1) "bin rwr" at first, (2) then  "hic interaction postprocess"
prefix="dataset_name"                                    #this will be used as a prefix for output file names
method="sliding_window"                                  #method for RWR computation, use inverse for chromosome-wide computation and "sliding_window" for faster computation

############################################################################


if [[ "$parallelism" == "parallel" ]]; then
	mpirun -np $number_of_processors python $snapHiC_dir/snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps --prefix $prefix --parallel
elif [[ "$parallelism" == "threaded" ]]; then
	python $snapHiC_dir/snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps  --prefix $prefix --threaded -n $number_of_processors
else
	python $snapHiC_dir/snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file  --prefix $prefix --steps $steps
fi
