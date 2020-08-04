<img src="nanoTRF.png" width="550" >

# NanoTRF: software tool to *de novo* search high-copy tandem repeats in Oxford Nanopore Technologies (ONT) plant DNA sequencing data



## Table of Contents

- [Introduction](#introduction)
- [Getting Started](#getting) 
  - [Building  nanoTRF from  source files](#building)
- [Commands and options](#cmd)
- [Input](#input_output)
- [Output](#output)
- [Running nanoTRF](#running)
  - [Usage](#usage)
- [Authors](#authors)
- [License](#license)
## <a name="getting"></a>Getting Started
### <a name="building"></a>Building  nanoTRF from  source files

Download the [latest release](https://github.com/Kirovez/nanoTRF/releases):
```
wget https:/https://github.com/Kirovez/nanoTRF/releases/download/v1.0.0/nanoTRF-v1.0.0.tar.gz
tar -zxvf nanoTRF-v1.0.0.tar.gz && cd TideHunter-v1.0.0

```
**Prerequisites**
nanoTRF requires:

- blastn and makeblastdb programs. The paths to these programs can be set via -bn and -mb flags, respectively
- TideHunter programm. It is recommended to download the [latest release of TideHunter](https://github.com/yangao07/TideHunter/releases).The paths to these programs can be set via **-pTH** flags
- Canu programm. The latest release [can be download here](http://github.com/marbl/canu/releases). The paths to these programs can be set via **-cu** flags
- python >= v3.6 python packages to be installed: biopython, networkx (run command: pip install matplotlib biopython networkx)


## <a name="introduction"></a>Introduction

NanoTRF is software tool to *de novo* search high-copy tandem repeats which is designed for raw long-read sequnces.

It works with Oxford Nanopore Technologies (ONT) sequencing data

### <a name="cmd"></a>Command and options

**-h, --help**  - show this help message and exit

**-r,--reads** - path to FastQ or Fasta file **(required argument!!!)**

**-out,--out_directory** - path to work directory for output files where will be saved **(required argument!!!)**

**-pTH, --path_TH** - path to the location of TideHunter **(required argument!!!)**

**-cu,--canu**  - path to the location of the Canu **(required argument!!!)**

**-bn,--blast**  - path to blastn executabled". ***Default='blastn'***

**-mb,--makedb**  - path to makeblastdb executable. ***Default='makeblastdb'***

**-w, --wordsize** - word size for wordfinder algorithm (length of best perfect match). ***Default = 22***

***-ev, --evalue*** -  expectation value (E) threshold for saving hits. ***Default = 2***

**-m,--max_abundancy**  - the proportion of amount lengths all tandem repeats in one cluster to length all the reads. ***Default = 0.0001***
                        
**-cons, --consensus_name** - file name with consensus sequences. ***Default='consensus.fasta'***

**-th, --threads**  - number of threads for running Blast. ***Default = 4***

**-lg, ---log_file**  - this file list analysis parameters, modules and files,contains messages generated 
on the various stages of the NanoTRF work. It allows tracking events that
happens when NanoTRF runs. Default - loging.log ***Default = 'loging.log'***

**-mOVe, --min_Overlap** - number of overlapping nucleotides between repeats in one cluster. ***Default = 15***

**-del, --opt_delete** - remove unnecessary large files and directories from working directory. ***Default='d'*** **d**-(delete files and directories), 
**c**-save all files and directories

## <a name="input_output"></a>Input
NanoTRF works with FASTA and FASTQ formats.

## <a name="output"></a>Output

NanoTRF generates consensus sequences in FASTA format.
### <a name="running"></a>Running nanoTRF

#### <a name="usage"></a>Usage

To generate consensus sequences in FASTA format file (with usage default optional arguments):
```
python3 ./nanoTRF.py -r test.fasta -pTH ./bin/TideHunter -cu ./bin/canu ./test/
```
To generate consensus sequences in FASTA format file, change number of theads that will be used and remove all unnecessary files and directories:
```
python3 ./nanoTRF.py -r test.fasta -pTH ./bin/TideHunter -cu ./bin/canu ./test/ -th 30 -del c
```
## <a name="authors"></a>Authors
**Ilya Kirov**

**Elizaveta Kolganova**

## <a name="license"></a>License
This project is licensed under the **MIT** License



