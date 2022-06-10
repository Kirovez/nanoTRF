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

```
## <a name="input_output"></a>Input
NanoTRF works with FASTA and FASTQ formats.

## <a name="output"></a>Output
### <a name="output"></a>Tabular file
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

```

### <a name="output"></a>html file with cluster pictures

## <a name="authors"></a>Authors
**Elizaveta Kolganova** [liza.colg@gmail.com](liza.colg@gmail.com)

**Ilya Kirov** [kirovez@gmail.com ](kirovez@gmail.com )

## <a name="ackn"></a>Acknowledgement
The project was financially supported by Russian Foundation for Basic Research (RFBR project № 17-00-00336)

## <a name="license"></a>License
This project is licensed under the [MIT](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt) License



