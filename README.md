
# TROP : Tandem Repeat Overlap Parser

make intro blabla Tandem Repeat Finder (TRF) 
TROP is R script (not best langage and should be recoded in perl or others
list existing tools (copy from the msc the list)

"In previous work [3], we observed reads made of abundant tandem repeat artifact not representative of the genome, which often led to the assembly of artifactual tandem contigs. To date, the tandem output at the flow cell and channel level, and across different library preparation has not been carefully investigated."
"Finally, we used TRF as a classic tool for tandem detection [11] with a custom script [12]) to join tandems whose coordinates overlap. Future studies desiring to characterize tandem content could test recent software development for comparison. These include nucleotide-based detection softwares, such as TideHunter [13], NCRF [14], NanoSTR [15], mTR [16], and signal-based softwares, such as DeepRepeat [17] and WarpSTR [18], all of which are potentially computationally much faster than TRF [11]. Their use to develop trimming tools represent a potential avenue of research to remove/mask artifactual tandems from raw reads prior to assembly. Indeed, this may be important as we spotted the integration of artifactual tandem repeats on one of our chromosome assembly (MIP2461, Table 2, see indices). We hypothesize that such issue may happen when sufficient reads share the same artifactual tandem sequence artifact on their edge (Fig. 5)."
"The TRF report was then parsed to join overlapping tandems and compute total tandem length per read as well as their location with a custom script written in R [46] named TROP (Tandem Repeat Overlap Parser [12])."


# Prerequisites
**1) Running TRF to produce a ```.dat``` file**

For the present tutorial, TRF was run with the following options on a fasta file of R9.4 nanopore reads (fasta file ```NA_all.fasta``` not provided)

```trf NA_all.fasta 2 5 7 80 10 50 2000```

Please visit TRF's github repository for further information on TRF installation and running options (https://github.com/Benson-Genomics-Lab/TRF)

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
