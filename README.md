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
**Place all files necessary in this folder on your computer (above is an example, adjust to the folder path on your computer)**

**2) Load necessary R packagae (install it prior)**

```
require(data.table)
```
**3) Load the script (after you download it from this repository and placed it your working directory)**

```
source("TROP.R")
```

Tandem Repeat Finder ran with the following options on file ```NA_all.fasta```

```trf NA_all.fasta 2 5 7 80 10 50 2000```
