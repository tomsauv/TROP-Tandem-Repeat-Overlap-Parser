
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

This file reports the total tandem length detected per read (one entry per read)
```
seqname trgroup trlength seqlength trprop
ffca2289-5365-42e2-9474-0443943c14c3 1 153 557 27.47
ffcb3273-02ff-463a-9748-87067e73c33b 0 0 213 0
ffd00425-4921-41c7-bad9-d283a8ddb9ec 0 0 2704 0
ffd2ebb1-f0bd-43ec-9f64-be4d31b88c49 0 0 1545 0
ffd3cc81-79e6-483e-91b9-9286d69a9ef1 2 123 4182 2.94
ffd3e82e-47ce-4daa-b79a-9ce77c05f679 1 56 244 22.95
ffd7e1ec-6e9e-4437-a20b-84e9ea4be872 1 261 2549 10.24
ffd907b5-aec6-4a16-86d0-5304e2f92dd7 5 623 1108 56.23
...
```

```2_tandem_locations_per_read.txt```

This file reports the location of each tandem region, thus multiple entries per read may occur.

The tandem location is reported for both of the 5' edge (positive values) and 3' edge (as negative values) for plotting purposes (SEE TUTO PDF/WIKI).   

```
seqname trgroup start_5prime stop_5prime trlength seqlength trprop start_3prime stop_3prime
ffca2289-5365-42e2-9474-0443943c14c3 1 322 475 153 557 27.47 -83 -236
ffcb3273-02ff-463a-9748-87067e73c33b 0 0 0 0 213 0 0 0
ffd00425-4921-41c7-bad9-d283a8ddb9ec 0 0 0 0 2704 0 0 0
ffd2ebb1-f0bd-43ec-9f64-be4d31b88c49 0 0 0 0 1545 0 0 0
ffd3cc81-79e6-483e-91b9-9286d69a9ef1 1 31 124 93 4182 2.22 -4059 -4152
ffd3cc81-79e6-483e-91b9-9286d69a9ef1 2 148 178 30 4182 0.72 -4005 -4035
ffd3e82e-47ce-4daa-b79a-9ce77c05f679 1 28 84 56 244 22.95 -161 -217
ffd7e1ec-6e9e-4437-a20b-84e9ea4be872 1 17 278 261 2549 10.24 -2272 -2533
ffd907b5-aec6-4a16-86d0-5304e2f92dd7 1 26 377 351 1108 31.68 -732 -1083
ffd907b5-aec6-4a16-86d0-5304e2f92dd7 2 582 627 45 1108 4.06 -482 -527
ffd907b5-aec6-4a16-86d0-5304e2f92dd7 3 694 722 28 1108 2.53 -387 -415
ffd907b5-aec6-4a16-86d0-5304e2f92dd7 4 817 973 156 1108 14.08 -136 -292
ffd907b5-aec6-4a16-86d0-5304e2f92dd7 5 1010 1053 43 1108 3.88 -56 -99
```

```3_tandem_output_summary.txt```

This file reports the number of sequences containing detected tandem and and overall tandem bp. Proportions are also computed:

e.g. below ~48% of sequences contained tandems and the detected tandem length representied >12% of total base pairs. 
```
class nseq nseq_perc bp bp_perc
no_tandem 12099 51.96 31589056 87.3
tandem 11185 48.04 4593376 12.7
total 23284 100 36182432 100
```

Fields meaning for the three files are as follows:

**seqname**: sequence name

**trgroup**: tandem region identifier in a sequence read/contig read (as increasing counts, e.g. region 1, 2, 3, etc)

**trlength**: length of the reported tandem repeat region 

**seqlength**: length of the sequence read/contig

**trprop**: tandem repeat region proportion as compared to the read/contig total length

**start_5prime**: starting location of the tandem repeat region measured from the 5' edge

**stop_5prime**: stopping location of the tandem repeat region measured from the 5' edge

**start_3prime**: starting location of the tandem repeat region measured from the 3' edge (reported as a negative value)

**stop_3prime**: stopping location of the tandem repeat region measured from the 3' edge (reported as a negative value)

**nseq**: number of sequence reads/contigs

**nseq_perc**: percentage of sequence reads/contigs

**bp**: total base pairs for non-tandem or tandem repeat region

**bp_perc**: percentage of total base pairs counted for non-tandem and tandem repeat region


# How to cite this script

Sauvage T, Cormic A, Passerini D. A comparison of Oxford nanopore library strategies for bacterial genomics. BMC Genomics. 2023 [doi:10.1186/s12864-023-09729-z](https://doi.org/10.1186/s12864-023-09729-z)

or directly

Sauvage T. TROP : Tandem Repeat Overlap Parser. https://github.com/tomsauv/TROP-Tandem-Repeat-Overlap-Parser. Accessed September XX, 2023. DOI XXXXXXXXXXXXXXXXXXX

# Additional references

Benson G. Tandem Repeats Finder: a program to analyze DNA sequences. Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573

