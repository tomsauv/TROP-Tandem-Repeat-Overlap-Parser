# TROP : Tandem Repeat Overlap Parser


The script requires the ```.dat``` file produced by TRF. 

Here, we use the example file below from a nanopore flongle run

```
NA_all.fasta.2.5.7.80.10.50.2000.dat
``` 

**1) Set your working directory**

```
setwd("C:/Users/tomsauv/my_R_analyses")
``` 
Place all files necessary in this folder on your computer (above is an example, adjust to the folder path on your computer)

**2) Load necessary R packagae (install it prior)**

```
require(data.table)
```
**3) Load the script (after you download it from this repository and placed it your working directory)**

```
source("TROP.R")
```

# Running Tandem Repeat Finder (TRF) 

Further explanations, options for running TRF can be found at the authors dedicated repository below:
https://github.com/Benson-Genomics-Lab/TRF



For the above tutorial and provided files, Tandem Repeat Finder ran with the following options on the fasta file ```NA_all.fasta```

```trf NA_all.fasta 2 5 7 80 10 50 2000```


# How to cite TROP

Sauvage T, Cormic A, Passerini D. A comparison of Oxford nanopore library strategies for bacterial genomics. BMC Genomics. 2023 xxxxxxxxxx

or directly

Sauvage T. TROP : Tandem Repeat Overlap Parser. https://github.com/tomsauv/TROP-Tandem-Repeat-Overlap-Parser. Accessed September XX, 2023. DOI XXXXXXXXXXXXXXXXXXX

# Additional references

Benson G. Tandem repeats finder: a program to analyze DNA sequences. Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573
