<img src="nanoTRF.png" width="550" >

# NanoTRF: a pipeline for *de novo* identification and sequence assembly of high-copy tandem repeats in raw Oxford Nanopore plant DNA sequencing data


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
git clone https://github.com/Kirovez/nanoTRF.git
cd nanoTRF
```
**Prerequisites**
nanoTRF requires:

- **blastn** and **makeblastdb** programs. The paths to these programs can be set via **-bn** and **-mb** flags, respectively
- **TideHunter** programm. It is recommended to download the [latest release of TideHunter](https://github.com/yangao07/TideHunter/releases).The path to these program **must** be set via **-pTH** flag
- **Canu** programm. The latest release [can be download here](http://github.com/marbl/canu/releases). The path to Canu **must** be set via **-cu** flag
- **python >= v3.6**
- python packages to be installed: **biopython**, **networkx** To install these packages run the following command

```
pip install matplotlib biopython networkx python-louvain
```

or

```
pip3 install matplotlib biopython networkx python-louvain
```

**Important note!** If you have `community` python module installed you need to delete it because it interferes with `python-louvain` module used by nanoTRF. Use this command to delete `community` module:
```
pip3 uninstall community

```

## <a name="introduction"></a>Introduction

NanoTRF is software tool to *de novo* search high-copy tandem repeats which is designed for raw long-read sequnces.

It works with Oxford Nanopore Technologies (ONT) sequencing data

## <a name="cmd"></a>Command and options

### Required arguments***

**-r,--reads** - path to FastQ or Fasta file

**-out,--out_directory** - path to work directory for output files where will be saved

**-pTH, --path_TH** - path to the location of TideHunter

**-cu,--canu**  - path to the location of the Canu

### Optional arguments
**-h, --help**  - show this help message and exit
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

Test nanoTRF with the test files (nanoTRF/test_seq/test_seq.fa) using 30 threads:
```
python3 ./nanoTRF.py -r ./test_seq/test_4th_Linum.fasta -pTH ../TideHunter-v1.4.2/bin/TideHunter -cu ../canu/Linux-amd64/bin/canu -out ./test/
```
To generate consensus sequences in FASTA format file, change number of theads that will be used and remove all unnecessary files and directories:
```
python3 ./nanoTRF.py -r ./test_seq/test_4th_Linum.fasta -pTH ../TideHunter-v1.4.2/bin/TideHunter -cu ../canu/Linux-amd64/bin/canu -out ./test/ -th 30 -del c
```
## <a name="authors"></a>Authors

**Elizaveta Kolganova**
**Ilya Kirov**
## Acknowledgement
The project was financially supported by Russian Foundation for Basic Research (RFBR project № 17-00-00336)

## <a name="license"></a>License
This project is licensed under the **MIT** License



