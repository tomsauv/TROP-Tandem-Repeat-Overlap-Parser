
# TROP : Tandem Repeat Overlap Parser

# Why this script?

```trop``` was thought about after artifactual sequences made entirely of tandem repeats were observed in nanopore sequencing (Sauvage et al. 2023). In seeking to characterize and quantify this sequencing noise, we used a classic tool known as Tandem Repeats Finder (```trf```, Benson 1999).

The problem is that ```trf``` reports multiple overlapping tandem boundaries for a given repetitive region, which is due to the imperfect nature of repetitive DNA (added to sequencing errors...) and the period of the detected repetitive motif (linked to ```trf``` parameters).

**Thus, ```trop``` parses the lowest and highest boundary of a repetitive region to report the location of tandem region along a read and overall tandem repeat output**

![Screenshot](fig/trop_scheme.png)

```trop``` can parse the output ```trf``` produces from any ```fasta``` file containing tandem. However, here, we are more focussed on long raw sequence reads (nanopore, pacbio) or even assembled contigs

# Prerequisites
**1) Running ```trf``` to produce a ```.dat``` file**

For the present tutorial, ```trf``` was run with the following options on a ```fasta``` file of R9.4 nanopore reads (```NA_all.fasta```, not provided)

```trf NA_all.fasta 2 5 7 80 10 50 2000```, which produces ```NA_all.fasta.2.5.7.80.10.50.2000.dat```

If you open this provided example file in a simple text editor, results are displayed sequences by sequences as shown below:

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

There are a number of ways this can be done, for instance, you could use the ```awk``` command below (one-liner) on linux platform from the raw nanopore reads in ```fasta``` (below as ```NA_all.fasta```, file not provided):

```
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' NA_all.fasta | cut -d ' ' -f 1 | paste - - | sed -r 's/^>//g' > NA_all_tabulated_lengths.txt
```

Note that here, we split the nanopore sequence very long titles to only keep the main identifier and we output a file named ```NA_all_tabulated_lengths.txt```. 

The resulting tabulated file format is as below. The first column is the sequence name, the second column is the sequence length in bp.

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
Place all files necessary in this working directory folder on your computer (above is an example, adjust to the folder path and name on your computer)

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

Thus, with the provided example files, we would type:

```
trop("NA_all.fasta.2.5.7.80.10.50.2000.dat","NA_all_tabulated_lengths.txt")
```

# File output

trop produces three output file:

```1_tandem_length_per_reads.txt```
```
seqname trgroup trlength seqlength trprop
000329bc-d1c5-41e0-b788-33edb403e025 0 0 526 0
00058de9-6fa2-4997-8dcc-f127fb4919fd 0 0 2372 0
0005c83e-0e71-43e0-b1fc-7e76513ba182 0 0 647 0
000aadcb-5592-4c2f-ac6f-58cc38db29b9 0 0 3465 0
000ab28b-db24-4ee5-83d2-10d8fcd7a5e4 1 250 484 51.65
000ceda9-b66a-430c-9d0e-213c598e62c2 1 297 1793 16.56
000f9e01-abbb-4e97-99e1-2b7ef8cd862c 0 0 1739 0
0010bd89-2875-4da9-ac68-630723a64513 2 519 688 75.44
...
```

```2_tandem_locations_per_read.txt```
```
seqname trgroup start_5prime stop_5prime trlength seqlength trprop start_3prime stop_3prime
000329bc-d1c5-41e0-b788-33edb403e025 0 0 0 0 526 0 0 0
00058de9-6fa2-4997-8dcc-f127fb4919fd 0 0 0 0 2372 0 0 0
0005c83e-0e71-43e0-b1fc-7e76513ba182 0 0 0 0 647 0 0 0
000aadcb-5592-4c2f-ac6f-58cc38db29b9 0 0 0 0 3465 0 0 0
000ab28b-db24-4ee5-83d2-10d8fcd7a5e4 1 192 442 250 484 51.65 -43 -293
000ceda9-b66a-430c-9d0e-213c598e62c2 1 1445 1742 297 1793 16.56 -52 -349
000f9e01-abbb-4e97-99e1-2b7ef8cd862c 0 0 0 0 1739 0 0 0
0010bd89-2875-4da9-ac68-630723a64513 1 29 375 346 688 50.29 -314 -660
0010bd89-2875-4da9-ac68-630723a64513 2 454 627 173 688 25.15 -62 -235
```

```3_tandem_output_summary.txt```
```
class nseq nseq_perc bp bp_perc
no_tandem 12099 51.96 31589056 87.3
tandem 11185 48.04 4593376 12.7
total 23284 100 36182432 100
```

# How to cite this script

Sauvage T, Cormic A, Passerini D. A comparison of Oxford nanopore library strategies for bacterial genomics. BMC Genomics. 2023 xxxxxxxxxx

or directly

Sauvage T. TROP : Tandem Repeat Overlap Parser. https://github.com/tomsauv/TROP-Tandem-Repeat-Overlap-Parser. Accessed September XX, 2023. DOI XXXXXXXXXXXXXXXXXXX

# Additional references

Benson G. Tandem Repeats Finder: a program to analyze DNA sequences. Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573

