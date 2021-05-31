## Applying alternative methods to impute single cell contact frequency 
By default, SnapHiC uses the random walk with restart (RWR) algorithm to impute 10Kb bin resolution contact frequency for each single cell. Users can apply alternative methods, such as Higashi (https://github.com/ma-compbio/Higashi), to impute scHi-C data, and use the script *normalize_imputed_matrices.py* to calculate the normalized contact frequency (i.e., Z-score normalization and outlier trimming proposed by SnapHiC).

In the current version (v0.1.0), SnapHiC parameters are optimized based on the RWR-imputed single cell contact frequency. If users apply alternative imputation methods, we recommend to test different SnapHiC parameters to obtain optimal loop calling results.  

Use the following command to run the script:  
```python  
python normalize_imputed_matrices.py -i <input_filename> -o <output_directory> -t <input_filetype> -b <binsize> -d <distance>
```  
In SnapHiC default, *binsize* and *distance* values are 1e4 and 2e6, respectively. The *input_filetype* can be either of "npz", "bedpe", or "hdf5".

The npz files are numpy files stored in the binary format.  

The bedpe files should contain contact frequency in the upper triangular of the scHi-C contact matrix, i.e. the second bin coordinate should be larger than the first bin coordinate. These bedpe files should include all intra-chromosomal bin pairs within the pre-specified 1D genomic distance range (2Mb by default), even if the imputed contact frequency is 0. The matrix size is determined by the maximal value of the second bin coordinate in the bedpe file.

The hdf files should follow the output format described by Higashi (https://github.com/ma-compbio/Higashi). 

If the files are in the npz or bedpe format, users need to provide the directory of the files with the *-i* argument of the script. Otherwise, users need to provide the directory followed by the hdf file name.
