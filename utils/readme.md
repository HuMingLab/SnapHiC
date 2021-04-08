## Normalizing values acquired by other imputation methods:
The script *normalize_imputed_matrices.py* is provided for computing the normalized values (z-score normalization and outlier trimming proposed by SnapHiC) for the users who have performed other imputation methods instead of the RWR method described in the SnapHiC.   
Use the following command to run the script:  
```python  
python normalize_imputed_matrices.py -i <input_filename> -o <output_directory> -t <input_filetype> -b <binsize> -d <distance>
```  
Default *distance* and *binsize* values used in SnapHiC are 1e4 and 2e6. The *input_filetype* can be either of "hdf5", "bedpe", or "npz".  
The npz files are numpy files stored in binary format.  
Bedpe files should be upper triangular, i.e. the second bin coordinate on each row should be larger than the first bin coordinate. They should also include only intra-chromosomal binpairs and the value for the very last binpair should be stored as well even if it is 0. (The matrix size is determined by the maximum second coordinate in the bedpe file).  
If the input is in hdf format, it should follow the output format described by higashi. Namely, it should contain one dataset called "coordinates" storing the x and y coordinates for all the cells. For each cell, there should be a separate dataset in the same file. These additional datasets are called "cell_x" where x is an integer number (cell_0, cell_1, etc.). Each of the "cell_x" datasets contains one vector of imputed values, that will be matched with the coordinates in the "coordinates" dataset.  
If the files are in npz or bedpe format, you should only provide the directory of the files with the *-i* argument of the script. Otherwise, you should provide the path followed by the hdf filename.


