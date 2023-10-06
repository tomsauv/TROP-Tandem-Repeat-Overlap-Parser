
# TROP : Tandem Repeat Overlap Parser

# Why this script

```trop``` was thought about after artifactual sequences made entirely of tandem repeats were observed in nanopore sequencing (Sauvage et al. 2023). In seeking to characterize and quantify this type of sequencing noise, we used a classic tool known as Tandem Repeats Finder (```trf```, Benson 1999).

The problem is that ```trf``` reports multiple overlapping tandem solution/boundaries for a given repetitive region, which is due to the imperfect nature of repetitive DNA and the period of the repetitive motif (and ```trf``` parameters).

Thus, ```trop``` parses the lowest and highest boundary of a repetitive region to better characterise overall tandem content along a read (output and location) 

```trop``` can parse any ```fasta``` file ```trf``` will accept. Here we ar emore interested in long sequence reads (nanopore, pacbio) or even assembled contigs

# Prerequisites
**1) Running ```trf``` to produce a ```.dat``` file**

For the present tutorial, ```trf``` was run with the following options on a ```fasta``` file of R9.4 nanopore reads (```NA_all.fasta```, not provided)

```trf NA_all.fasta 2 5 7 80 10 50 2000```, which produces ```NA_all.fasta.2.5.7.80.10.50.2000.dat``` for later parsing

If you open this provided example file in a simple text editor, results sequences by sequences look as below:

```
Tandem Repeats Finder Program written by:

Gary Benson
Program in Bioinformatics
Boston University
Version 4.09


Sequence: e387a9f8-ebed-4d86-81a1-087c17c96745 runid=e3cf134550fca3678f03a505c0ce5509328f672c sampleid=no_sample read=36148 ch=109 start_time=2022-04-05T03:45:24Z model_version_id=2021-05-17_dna_r9.4.1_minion_



Parameters: 2 5 7 80 10 50 2000


59 185 52 2.4 54 79 11 168 28 14 25 30 1.95 CTTTCTGTTGGTGCTGATATTGCTGAAAATAGAGCGACAGGCAAGAGACAATAT CTTTCTGTTGGTGCTGATGGCTTGAAAATAGAGCGACAGGCAAGAGACAATATCTTTCATTGGTATGATATTGCTGAAGATAGAGCGACAGGCAAGACAATATCTTTCTGTTGGTGCAGATATTGCT


Sequence: d80c3c51-9115-419c-a8cf-3245e320a981 runid=e3cf134550fca3678f03a505c0ce5509328f672c sampleid=no_sample read=11233 ch=120 start_time=2022-04-04T14:00:13Z model_version_id=2021-05-17_dna_r9.4.1_minion_

...
```



Please visit ```trf```'s github repository for further information on ```trf``` installation, running and output options (https://github.com/Benson-Genomics-Lab/TRF)

**2) Preparing a tabulated sequence length file**

There are a number of ways this can be done, for instance, you could use the ```awk``` command below (one-liner) on linux platform from the raw nanopore reads in ```fasta```:

```
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' NA_all.fasta | cut -d ' ' -f 1 | paste - - | sed -r 's/^>//g' > NA_all_tabulated_lengths.txt
```

Note that here, we split the nanopore sequence title to only keep the main identifier. 

The resulting tabulated file (named above ```NA_all_tabulated_lengths.txt```) format is as below.

The first column is the sequence name, the second column is the sequence length in bp.

```
80800063-c64b-4493-a80c-6e2c6bfbe01a	823
a9fc5d29-acd8-4858-81c0-5cf6434708c0	894
66c1ef31-ac67-447a-ad27-4fdb74f47c74	2845
9cbcb1ea-5ae0-4333-a6f6-fcca75430742	3096
12ede4c6-f934-4566-bdd3-dc0e8efa004b	4402
3fa38391-78fc-47de-927d-1b806dc12ac9	456
7a5fe376-fa8b-4552-8803-04015cab89bd	2961
...
```

# Running TROP

**1) Open R and set your working directory**
```
setwd("C:/Users/tomsauv/my_R_analyses")
``` 
Place all files necessary in this folder on your computer (above is an example, adjust to the folder path on your computer)

**2) Load the ```data.table```  R package (install it prior)**
```
require(data.table)
```

**3) Load the script (after you download it from this repository and place it your working directory)**
```
source("trop.r")
```

**4) Run the script by filling the name of files to be parsed as below**

The first argument is the ```.dat``` file.

The second argument is the name of the sequence length file.


```
trop("NA_all.fasta.2.5.7.80.10.50.2000.dat","NA_all_tabulated_lengths.txt")
```

**OUTPUT**


# How to cite TROP

Sauvage T, Cormic A, Passerini D. A comparison of Oxford nanopore library strategies for bacterial genomics. BMC Genomics. 2023 xxxxxxxxxx

or directly

Sauvage T. TROP : Tandem Repeat Overlap Parser. https://github.com/tomsauv/TROP-Tandem-Repeat-Overlap-Parser. Accessed September XX, 2023. DOI XXXXXXXXXXXXXXXXXXX

# Additional references

Benson G. Tandem Repeats Finder: a program to analyze DNA sequences. Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573

# Comment

speed etc....perl
from a nanopore flongle run, small file provided
memory, made to use less memory
