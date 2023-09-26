
# TROP : Tandem Repeat Overlap Parser

make intro blabla Tandem Repeat Finder (TRF) 
TROP is R script (not best langage and should be recoded in perl or others
list existing tools (copy from the msc the list)

"In previous work [3], we observed reads made of abundant tandem repeat artifact not representative of the genome, which often led to the assembly of artifactual tandem contigs. To date, the tandem output at the flow cell and channel level, and across different library preparation has not been carefully investigated."
"Finally, we used TRF as a classic tool for tandem detection [11] with a custom script [12]) to join tandems whose coordinates overlap. Future studies desiring to characterize tandem content could test recent software development for comparison. These include nucleotide-based detection softwares, such as TideHunter [13], NCRF [14], NanoSTR [15], mTR [16], and signal-based softwares, such as DeepRepeat [17] and WarpSTR [18], all of which are potentially computationally much faster than TRF [11]. Their use to develop trimming tools represent a potential avenue of research to remove/mask artifactual tandems from raw reads prior to assembly. Indeed, this may be important as we spotted the integration of artifactual tandem repeats on one of our chromosome assembly (MIP2461, Table 2, see indices). We hypothesize that such issue may happen when sufficient reads share the same artifactual tandem sequence artifact on their edge (Fig. 5)."
"The TRF report was then parsed to join overlapping tandems and compute total tandem length per read as well as their location with a custom script written in R [46] named TROP (Tandem Repeat Overlap Parser [12])."

put graph? fig?

# Prerequisites
**1) Running TRF to produce a ```.dat``` file**

For the present tutorial, TRF was run with the following options on a fasta file of R9.4 nanopore reads (fasta file ```NA_all.fasta``` not provided)

```trf NA_all.fasta 2 5 7 80 10 50 2000```, which produces ```NA_all.fasta.2.5.7.80.10.50.2000.dat``` for later parsing by TROP

Please visit TRF's github repository for further information on TRF installation and running options (https://github.com/Benson-Genomics-Lab/TRF)

**2) Preparing a tabulated sequence length file**

There are a number of ways this can be done, for instance, you could use the awk command below (one-liner) on linux platform with the raw fasta reads:

```
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' NA_all.fasta | cut -d ' ' -f 1 | paste - - | sed -r 's/^>//g' > NA_all_tabulated_lengths.txt
```

Note that here, we split the nanopore sequence title to only keep the main identifier. The resulting tabulated file should be as below to be processed by TROP
(The first column is the sequence name, second column is the sequence length in bp)

```
80800063-c64b-4493-a80c-6e2c6bfbe01a	823
a9fc5d29-acd8-4858-81c0-5cf6434708c0	894
66c1ef31-ac67-447a-ad27-4fdb74f47c74	2845
9cbcb1ea-5ae0-4333-a6f6-fcca75430742	3096
12ede4c6-f934-4566-bdd3-dc0e8efa004b	4402
3fa38391-78fc-47de-927d-1b806dc12ac9	456
7a5fe376-fa8b-4552-8803-04015cab89bd	2961
eb57c2a6-e21d-4eb4-bdcd-4ecd60c74728	1254
42e88fd2-8089-411f-a8e2-a58a2c36c8e9	1288
86f082f7-68a7-4c91-871c-d3d9608e231f	719
8e6f5c79-6e2e-40cc-b785-bcf6a7316b4f	2074
239d894a-0adf-43a2-9c76-e6b68910b21a	3940
...
```





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
