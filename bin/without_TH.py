import logging
import os
from Bio import SeqIO
from bin.helpers.help_functions import getLog
class without_TH():
    def __init__(self, outTab, outFasta, log_file):

        self.outFasta = outFasta
        self.outTab = outTab
        self.log_th = getLog(log_file, "without TH")
        self.tab2fasta()
    def tab2fasta(self):
        cnt = 0
        with open(self.outFasta, 'w') as outFasta, open(self.outTab) as outTab:
            for lines in outTab:
                sp = lines.rstrip().split("\t")
                seq = ">{0}*{1}*{2}*{3}\n{4}\n".format(sp[0],sp[1], sp[5], sp[6], sp[10])
                outFasta.write(seq)
                cnt += 1
        logging.info("Number of tandem repeats is: {}".format(cnt))
