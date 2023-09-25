
# TROP : Tandem Repeat Overlap Parser

make intro blabla Tandem Repeat Finder (TRF) 
TROP is R script (not best langage and should be recoded in perl or others
list existing tools (copy from the msc the list)

# Prerequisites
**1) Running TRF to prodcue a ```.dat``` file**

For the present tutorial and provided files, TRF ran with the following options on the (provided) fasta file ```NA_all.fasta```
(nanopore reads)

```trf NA_all.fasta 2 5 7 80 10 50 2000```

Further explanations, options for running TRF can be found at the authors dedicated github repository below:
https://github.com/Benson-Genomics-Lab/TRF

**2) Preparing a tabulated sequence length file**

how to do etc...for nanopore

# Running TROP

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
**3) Load the above ```.dat``` file with R function ```scan```**
```
mydf <- scan("NA_all.fasta.2.5.7.80.10.50.2000.dat", what = "", sep = "\n")
```

**3) Load the script (after you download it from this repository and place it your working directory)**
```
source("TROP.R")
```




# How to cite TROP

Sauvage T, Cormic A, Passerini D. A comparison of Oxford nanopore library strategies for bacterial genomics. BMC Genomics. 2023 xxxxxxxxxx

or directly

Sauvage T. TROP : Tandem Repeat Overlap Parser. https://github.com/tomsauv/TROP-Tandem-Repeat-Overlap-Parser. Accessed September XX, 2023. DOI XXXXXXXXXXXXXXXXXXX

# Additional references

Benson G. Tandem Repeats Finder: a program to analyze DNA sequences. Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573
