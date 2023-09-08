# weighted_average_stdev.py #
Juan F. Santana, Ph.D. (<juan-santana@uiowa.edu>), University of Iowa, Iowa City, I.A.

This software computes the weighted mean location, standard deviation, and count of the 5', 3' termini or fragment centers within a specified region. It accepts multiple replicates as input, each with its associated correction values. The reads analyzed can be strand specific if specified.

# File requirements #
The input regions file should be a six-column, tab-delimited bed file containing chromosome, start and end positions, and the strand information for each region. The regions can be of any length as long as it is an even number and the center is a feature under study (e.g. transcription start site). 
 
| chr6 | 142946246 | 142946446 | Gene_A | 255 | - |
|:----:|:---------:|:---------:|:------:|:---:|:-:|

The input fragments file should be a six-column, tab-delimited bed file containing chromosome, start and end positions, and the strand information for each fragment.

| chr6 | 142946247 | 142946298 | A00876:119:HW5F5DRXX:2:2207:29170:1157 | 255 | - |
|:----:|:---------:|:---------:|:--------------------------------------:|:---:|:-:|


# Behavior #

This software serves as a computational tool designed to analyze and process fragment data within user-defined regions.

** Fragment Analysis**: Inspects the 5', 3' termini or the centers of fragments.
** Computational Operations**: Determination of the weighted mean location and standard deviation of these fragments.
** Count Analysis**: Quantifies the number of such fragments.
** Multiple Data Set Handling**: Processes multiple replicates simultaneously, adjusting calculations based on provided correction values for each replicate.
** User-Driven Boundaries**: Operates within user-specified boundaries.


# Dependencies #
### Python libraries ###
Pandas: https://pypi.org/project/pandas/

### External program ###
bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html, developed by the Quinlan laboratory at the University of Utah. 

# Example command usage #
```
python3 weighted_average_stdev.py plusminus1000_from_TSS_1000genes.bed \
                  -f PolII-DFF-ChIP-Seq-Rep1.bed PolII-DFF-ChIP-Seq-Rep2.bed PolII-DFF-ChIP-Seq-Rep3.bed \
                  -c 1.2 1.5 \
                  -o /home/user/dir/ \
                  -r -100 250 \
                  -t 5 \
                  -s \
```
# Parameter description #
```
regions: <str> Bed file of genomic regions of chosen length with the format described above

-f: <str> Bed file(s) of fragment positions with the format described above. Multiple data files (replicas) can be added at the same time

-c: <int | float> Correction factors for each data file in `-f` in the same order

-r: <int> <int> Range relative to the TSS to analyze, for exmaple -r -100 250

-t: <str> Choose from '5' or '3' (for transcriptional data) or 'centers' (for ChIP-seq data). The input should be full fragments, the program will automatically calculate the 5 or 3 prime ends or centers

-s: If invoked, stranded data will be used

-o: <str> Path to output file

```

