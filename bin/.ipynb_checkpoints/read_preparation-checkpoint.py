from Bio import SeqIO
import logging
from bin.helpers.Errors import InputError
#logging.basicConfig(filename='example.log', level=logging.DEBUG)

class PrepareReads():
    def __init__(self, read_file):
        logging.debug("Read preparation....")
        self.read_file = read_file
        self.cnt = 0
        self.total_bases = 0
        self.format = 'fasta'
        self.prepareReads()

    def _checkReadFormat(self, read_file):
        with open(read_file) as reads:
            first_line = reads.readline()
            if first_line.startswith('@'):
                logging.info('Input file is fastq')
                return 'fastq'
            elif first_line.startswith(">"):
                logging.info('Input file is fasta')
                return 'fasta'
            else:
                raise InputError("The format (fastq or fasta) of input file can not be determined! \n"
                                 "Ensure that first line starts with '>' (for fasta) or '@' (for fastq)")

    def _fq2fasta(self):
        outFasta = self.read_file + '.fasta'
        with open(outFasta,'w') as outFile:
            for seq in SeqIO.parse(self.read_file, "fastq"):
                SeqIO.write(seq, outFile, 'fasta')
                self.cnt += 1
                self.total_bases += len(seq.seq)
        print("Number of reads: {0}".format(self.cnt))
        logging.info("Input file was converted into fasta\n Number of reads in the input file: {0}".format(self.cnt))
        self.read_file = outFasta


    def _countTotal(self):
        for seq in SeqIO.parse(self.read_file, self.format):
            self.cnt += 1
            self.total_bases += len(seq.seq)
        print("Number of reads: {0}".format(self.cnt))
        logging.info("Number of reads in the input file: {0}".format(self.cnt))


    def prepareReads(self):
        format = self._checkReadFormat(self.read_file)
        if format == 'fastq':
            self.format = format
            self._fq2fasta()
        else:
            self._countTotal()

    # count number of sequences and
#PrepareReads(r'C:\Users\ikirov\Google Диск\Bioinformatics_2018\practice\seq.fastq')