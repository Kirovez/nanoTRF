<img src="nanoTRF.png" width="550" >

# NanoTRF: a software tool to *de novo* search high-copy tandem repeats in Oxford Nanopore Technologies (ONT) plant DNA sequencing data

Download the [latest release](https://github.com/Kirovez/nanoTRF/releases):
```
wget https:/https://github.com/Kirovez/nanoTRF/releases/download/v1.0.0/nanoTRF-v1.0.1.tar.gz
tar -zxvf nanoTRF-v1.0.0.tar.gz && cd TideHunter-v1.0.0
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

NanoTRF is a software tool to *de novo* search high-copy tandem repeats which is designed for raw long-read sequences. It works with Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) sequencing data.

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

Before you start, you need to make sure that all program and packages specified below is already installed on your computer. For running nanoTRF  you will need to specify the path of the program through special flags:

- blastn and makeblastdb programs. The paths to these programs can be set via **-b** and **-mb** flags, respectively
- TideHunter program. It is recommended to download the [latest release of TideHunter](https://github.com/yangao07/TideHunter/releases). The paths to these programs can be set via **-pTH** flags
- Canu program. The latest release [can be download here](http://github.com/marbl/canu/releases). The paths to these programs can be set via **-cu** flags
- java
- python >= v3.6
- python packages to be installed: biopython, networkx. To install these packages run the following command
```
 pip install matplotlib biopython networkx python-louvain
```
or
```
pip3 install matplotlib biopython networkx python-louvain
```
Important note! If you have a community python module installed you need to delete it because it interferes with the python-louvain module used by nanoTRF. Use this command to delete the community module:
```
pip3 uninstall community
```


## <a name="usage"></a>Usage

To generate consensus sequences in FASTA format file (with usage default optional arguments):
```
python3 ./nanoTRF.py -r test.fasta -pTH  -cu ./bin/canu -o./test/
```
To generate consensus sequences in FASTA format file, change the number of threads that will be used and remove 
all unnecessary files and directories (with usage TideHunter files) using 30 threads:
```
python3 ./nanoTRF.py -r test.fasta --cu ./bin/canu -o ./test/ -th 30 -d -T TH.tab TH.out.fasta
```
## <a name="cmd"></a>Command and options
```

Options:
  General options:
      -h --help               show this help message and exit
 

  Input:
    -r --reads          STR      path to FastQ or Fasta file (required argument!!!)
    -T --run_th         STR      path to output files of the TideHunter (if previously TideHunter was running by user): 
                                 table file with consensus sequences and fasta file with unique tandem repeats
  Scoring parameters for partial order alignment:
    -w --wordsize       INT      word size for wordfinder algorithm (length of best perfect match) (Default = 22)
    -w_f --wordsize_f   INT      word size for wordfinder algorithm (length of best perfect match) in 
                                 the Reclusting module (Default=15)
    -ev --evalue        INT      expectation value (E) threshold for saving hits (Default=2)

  Clustering parameters:
    -m --min_copy       INT      The minimum number of TRs copy in the data (Default=100)
    -mOVe --min_Overlap STR      the number of overlapping nucleotides between repeats in one cluster (De10]
    -ca --perc_abund    STR      minimum value of the TR cluster abundancy. ***Default = 0.009***

  Path to programm for running nanoTRF:
    -pTH --path_TH      STR      path to the location of TideHunter [TideHinter]
    -cu --canu          STR      path to the location of Canu (required argument!!!It's missing in the conda)
    -trf --TRF_run      STR      path to the location [trf]
    -b --blast          STR      path to blastn executabled [blastn]
    -mb --makedb        STR      path to makeblastdb executable [makeblastdb]

  Output:
    -o --out_directory  STR      path to work directory for output files where will be saved **(required argument!!!)
    -lg --log_filepath  STR      path to file which list analysis parameters, modules, and files, contains messages generated 
                                 in the various stages of the work [loging.log]
    -nano --nano_trf    STR      fasta file with the TRs consensus sequences [nanoTRF.fasta]
    -tab --nano_tab     STR      table file with the TRs abundancy [TR_info.tab]

  Сomputational resources:
    -th, --threads      STR      number of threads for running blast, canu. [4]

  Additional option:
    -d --dir_cleanup    STR      remove unncessary large files and directories from working directory [False]
    
    
-h, --help  - show this help message and exit

```
## <a name="input_output"></a>Input
NanoTRF works with FASTA and FASTQ formats.

## <a name="output"></a>Output
### <a name="output"></a>Tabular file
NanoTRF generates output in tabular format:
| №   | Column name | Description | 
|:---:|   :---      | ---        |
|  1  | Cluster     | Name and cluster number |
|  2  | TRs length  | Length of the TRs consensus sequence |
|  3  | Abundance   | Cluster abundancy (%)



### <a name="output"></a>Fasta file

NanoTRF generates TRs consensus sequences in FASTA format which contents information about TRs. The sequence descriptions have the following format:
```
>clustname monomer_length cluster_abund

clustname          cluster number (for example clust0)
monomer_length     length of the TRs sequence
cluster_abund      cluster abundancy
```
## <a name="authors"></a>Authors
**Elizaveta Kolganova** [liza.colg@gmail.com](liza.colg@gmail.com)

**Ilya Kirov** [kirovez@gmail.com ](kirovez@gmail.com )

## <a name="ackn"></a>Acknowledgement
The project was financially supported by Russian Foundation for Basic Research (RFBR project № 17-00-00336)

## <a name="license"></a>License
This project is licensed under the [MIT](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt) License



