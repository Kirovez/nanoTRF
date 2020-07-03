"""
This script takes path with fasta file, name of the oufile and number of threshold value that is required to add a particular atom  for searching consensus.
 Results file are fasta file with aligned sequences and file with consensus.
Example:
    /path/lk/ 0.5 outfile.txt
"""

from Bio.Align.Applications import MafftCommandline
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
import os
import argparse
import random


# function for running mafft - a multiple sequence alignment program
def maffRun(dir_fasta):
    # create a specific directory. if there exists, name of subsequent directory will be change through counter
    name_dir = '{1}{0}new_Dir/'.format(random.randint(0, 10000), dir_fasta)
    if os.path.exists(name_dir):
        create_dir = os.mkdir(name_dir)
    else:
        # create a specific directory for fasta file after mafft
        create_dir = os.mkdir(name_dir)

        for fasta in os.listdir(dir_fasta):
            # take only fasta files in directory
            if 'fasta' in fasta:
                # running mafft
                mafft_cline = MafftCommandline(input=fasta)
                stdout, stderr = mafft_cline()
                # name for fasta file with aligned sequences
                nameAl = '{0}alm.fasta'.format(fasta)
                FileFullPath = os.path.join(name_dir, nameAl)
                # write to the file of specific directory aligned sequnces
                with open(FileFullPath, 'w') as handle:
                    handle.write(stdout)
    return name_dir


# function for searching consensus for aligned sequences and write in file in fasta format
def alignIf(dir_fasta, value_trd):
    consensus_list = []
    for fasta in os.listdir(dir_fasta):
        print(fasta)
        if 'fasta' in fasta:
            FileFullPath = os.path.join(dir_fasta, fasta)
            alignment = AlignIO.read(open(FileFullPath), 'fasta')
            summary_align = AlignInfo.SummaryInfo(alignment)
            # definition of the consensus
            consensus = summary_align.dumb_consensus(threshold=float(value_trd), require_multiple=1)
            record = '>{0}'.format(fasta.split('.fasta')[0])
            consensus_list.append('>{0}{1}{2}'.format(record, '/n', consensus))

    return consensus_list


def writeFile(end_table, outFile, dir_fasta):
    FileFullPath = os.path.join(dir_fasta, outFile)
    with open(FileFullPath, 'w') as outfile:
        for key in end_table:
            outfile.write(str(key) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Searching consensus for different TR in clusters')
    parser.add_argument('directory', help="Directory with file")
    parser.add_argument('outfile', help="Outfile fasta with consensus")
    parser.add_argument('value_trd', help="Value of threshold for alignment")

    args = parser.parse_args()
    mafft_file = maffRun(args.directory)
    consensus_list = alignIf(mafft_file, args.value_trd)
    writeFile(consensus_list, args.outfile, mafft_file)
