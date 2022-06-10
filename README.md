<img src="nanoTRF.png" width="550" >

# NanoTRF: a software tool to *de novo* search high-copy tandem repeats in Oxford Nanopore Technologies (ONT) plant DNA sequencing data

Download the [latest release](https://github.com/Kirovez/nanoTRF/releases):
```
wget https://github.com/Kirovez/nanoTRF/archive/refs/tags/nanoTRF.tar.gz
tar -zxvf nanoTRF.tar.gz
```
 

## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing nanoTRF via conda](#conda)
  - [Building nanoTRF from source files](#building)
- [Getting Started](#getting) 
- [Usage](#usage)
- [Commands and options](#cmd)
- [Input](#input_output)
- [Output](#output)
- [Authors](#authors)
- [Acknowledgement](#ackn)
- [License](#license)

## <a name="introduction"></a>Introduction

NanoTRF is a software tool to *de novo* search high-copy tandem repeats designed for raw long-read sequences. It works with Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) sequencing data.

## <a name="install"></a>Installation

### <a name="conda"></a>Installing nanoTRF via conda
On Linux/Unix, nanoTRF can be installed via creating an environment from an environment.yml file:
 ```
 conda env create -f nanoTRF.yml
 ```
For running nanoTRF, please activate the conda environment:
 ```
 conda activate nanoTRF
 ```
 Your environment is ready to be used!
 
### <a name="conda"></a>Pre-built binary executable file for Linux/Unix

If you meet any issue with creating the environment, please try the pre-built binary file:

```
wget https:/https://github.com/Kirovez/nanoTRF/releases/download/v1.0.0/nanoTRF-v1.0.0.tar.gz
tar -zxvf nanoTRF-v1.0.0.tar.gz && cd TideHunter-v1.0.0
```
**We recommended create a folder and running the pipeline in the previously established directory**

Before you start, you need to ensure that all programs and packages specified below are already installed on your computer. For running nanoTRF, you will need to select the path of the program through special flags:

- blastn and makeblastdb programs. Users can set the paths to these programs via **-b** and **-mb** flags, respectively
- TideHunter program. It is recommended to download the [latest release of TideHunter](https://github.com/yangao07/TideHunter/releases). Users can set the paths to these programs via **-pTH** flags
- Canu program. The latest release [can be download here](http://github.com/marbl/canu/releases). The paths to these programs can be set via **-cu** flags
- java
- python >= v3.6
- Python packages to be installed: biopython, networkx. To install these packages, run the following command
```
 pip install matplotlib biopython networkx python-louvain
```
or
```
pip3 install matplotlib biopython networkx python-louvain
```
Important note! Suppose you have a community python module installed. In that case, you need to delete it because it interferes with the python-louvain module used by nanoTRF. Use this command to delete the community module:
```
pip3 uninstall community
```


## <a name="usage"></a>Usage

To generate consensus sequences in FASTA format file (with usage default optional arguments):
```
python3 ./nanoTRF.py -r test.fasta  -o./test/
```
If TideHunter output table (run with -f option) was generated before then you can pass this file via -T option and nanoTRF will skip TideHunter step

```
python3 ./nanoTRF.py -r test.fasta -o ./test/ -T TH.tab
```

## <a name="cmd"></a>Command and options
```

usage: nanoTRF.py [-h] [-r READS] [-pTH PATH_TH] [-T RUN_TH] [-cap CAP3] [-diamond DIAMOND] [-o OUT_DIRECTORY]
                  [-b BLAST] [-mb MAKEDB] [-w WORDSIZE] [-w_f WORDSIZE_F] [-ev EVALUE] [-mid MIN_ID]
                  [-bld QUERY_SBJ_LENGTH_DIFFERENCES_ALLOWED] [-mad MIN_ABUNDANCY_TO_DRAW] [-m MIN_COPY]
                  [-nano NANO_TRF] [-tab NANO_TAB] [-rexdb_fasta REXDB_FASTA] [-rexdb_tab REXDB_TAB] [-th THREADS]
                  [-lg LOG_FILE] [-mOVe MIN_OVERLAP] [-ca PERC_ABUND] [-c] [-maskws MASK_BLAST_WORD_SIZE]
                  [-maskcov MASK_BLAST_QUERY_COVERAGE] [-maskiden MASK_BLAST_IDENTITY]

A tool to clustering sequences in fasta file and searching consensus among the many sequences for each cluster

optional arguments:
  -h, --help            show this help message and exit
  -r READS, --reads READS
                        Path to FastQ or Fasta file
  -pTH PATH_TH, --path_TH PATH_TH
                        Path to the location of the TideHunter
  -T RUN_TH, --run_th RUN_TH
                        If you do not want to run TideHunter again and you have table file (-f 2 option in Tide
                        Hunter), type the path to this file here
  -cap CAP3, --cap3 CAP3
                        Path to the location of the Cap3
  -diamond DIAMOND, --diamond DIAMOND
                        Path to the location of DIAMOND
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Path to work directory for output files where will be saved
  -b BLAST, --blast BLAST
                        Path to blastn executabled
  -mb MAKEDB, --makedb MAKEDB
                        Path to makeblastdb executable
  -w WORDSIZE, --wordsize WORDSIZE
                        Word size for wordfinder algorithm (length of best perfect match)
  -w_f WORDSIZE_F, --wordsize_f WORDSIZE_F
                        Word size for Reblusting(length of best perfect match)
  -ev EVALUE, --evalue EVALUE
                        Expectation value (E) threshold for saving hits
  -mid MIN_ID, --min_id MIN_ID
                        minimum identity between monomers to be selected for clustering
  -bld QUERY_SBJ_LENGTH_DIFFERENCES_ALLOWED, --query_sbj_length_differences_allowed QUERY_SBJ_LENGTH_DIFFERENCES_ALLOWED
                        maximum differences in length between query and subject
  -mad MIN_ABUNDANCY_TO_DRAW, --min_abundancy_to_draw MIN_ABUNDANCY_TO_DRAW
                        Minimum genome abundancy for cluster of repeats to be drawn
  -m MIN_COPY, --min_copy MIN_COPY
                        The minimum number of TRs copy in the data
  -nano NANO_TRF, --nano_trf NANO_TRF
                        File name with consensus sequences, default name - nanoTRF.fasta
  -tab NANO_TAB, --nano_tab NANO_TAB
                        Table file with the TRs abundancy
  -rexdb_fasta REXDB_FASTA, --rexdb_fasta REXDB_FASTA
                        Fasta file with the RExDB protein sequences
  -rexdb_tab REXDB_TAB, --rexdb_tab REXDB_TAB
                        Table file with the RExDB annotation
  -th THREADS, --threads THREADS
                        Number of threads for running the module Blast and TideHunter
  -lg LOG_FILE, ---log_file LOG_FILE
                        This file list analysis parameters, modules and files, contains messages generated on the
                        various stages of the NanoTRF work. It allows tracking events that happens when NanoTRF runs.
                        Default =loging.log
  -mOVe MIN_OVERLAP, --min_Overlap MIN_OVERLAP
                        Number of overlapping nucleotides between repeats in one cluster
  -ca PERC_ABUND, --perc_abund PERC_ABUND
                        Minimum value of the TR cluster abundancy
  -c, --cleanup         Remove unncessary large files and directories from working directory
  -maskws MASK_BLAST_WORD_SIZE, --mask_blast_word_size MASK_BLAST_WORD_SIZE
                        word size of blastn masking of raw reads by cluster contig sequences
  -maskcov MASK_BLAST_QUERY_COVERAGE, --mask_blast_query_coverage MASK_BLAST_QUERY_COVERAGE
                        query (contig sequence) coverage in blastn masking of raw reads by cluster contig sequences
  -maskiden MASK_BLAST_IDENTITY, --mask_blast_identity MASK_BLAST_IDENTITY
                        minimum identity between query (contig sequence) and raw reads in blastn masking of raw reads
                        by cluster contig sequences
```
## <a name="input_output"></a>Input

NanoTRF works with FASTA and FASTQ formats.

## <a name="output"></a>Output

### <a name="output"></a>Tabular file `clust_abund.tab`

NanoTRF generates output in tabular format:
| №   | Column name | Description | 
|:---:|   :---      | ---        |
|  1  | Cluster     | Name and cluster number |
|  2  | min.Contig.Cap3.Length  | Min length of the contigs assembled by Cap3 |
|  3  | max.Contig.Cap3.Length   | Max length of the contigs assembled by Cap3
|  4  | Genome.portion   | Cluster abundancy in the genome (%)
|  5  | Contig1.sequence   | Sequence of Contig1 from consensus.fasta
|  6  | Subrepeats.seq   | Sequences of any detected subrepeats in Contig 1 by second run of TideHunter 
|  7  | Subrepeats.len   | Length of subrepeats sequences
|  8  | Annotation  | Transposon domains and number of reads with similarity to


### <a name="output"></a>Fasta file

NanoTRF generates 'consensus.fasta' file which contains TRs consensus sequences assembled by Cap3. 


### <a name="output"></a>Html file 

This file containes the information from tabular file and some pictures including graph layout, read coverage histogram and read coverage pie chart

## <a name="authors"></a>Authors

**Ilya Kirov** [kirovez@gmail.com ](kirovez@gmail.com )

**Elizaveta Kolganova** [liza.colg@gmail.com](liza.colg@gmail.com)



## <a name="ackn"></a>Acknowledgement
The project was financially supported by Russian Foundation for Basic Research (RFBR project № 17-00-00336)

## <a name="license"></a>License
This project is licensed under the [MIT](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt) License



