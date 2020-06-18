"""
makeblastdb -dbtype nucl -in nanoTRF_0.1M.fasta -out nanoTRF_0.1M.fasta

blastn -query merged_TR_rank_all.fasta -outfmt 6 -db  nanoTRF_0.1M.fasta -out merged_TR_rank_all_vs_nanoTRF_0.1M.out -window_size 22 -num_threads 100 -evalue 10
"""

from bin.helpers.help_functions import getLog
import os

class run_BLAST():
    def __init__(self, inFile, outFile, threads, wordsize, evalue, log_file):
        self.inFile, self.outFile, self.threads, self.wordsize, self.evalue = inFile, outFile, threads, wordsize, evalue
        self.bl_log = getLog(log_file, "BLAST module")
        self.edge_list_file = outFile + "edges.list"
        self.main()

    def filterOut_table(self):
        """
        :return: list of edges [(query, hit),(),...]
        """
        edge_cnt = 0
        with open(self.outFile) as inFile, open(self.edge_list_file, 'w') as outEdgeList:
            for lines in inFile:
                sp = lines.split("\t")
                if sp[0] != sp[1] and float(sp[10]) < 0.00001:
                    outEdgeList.write("{0}\t{1}\t{2}\n".format(sp[0], sp[1], sp[10]))
                    edge_cnt += 1

        print("NUmber of edges", edge_cnt)
        self.bl_log.info("NUmber of edges: {}".format(edge_cnt))


    def main(self):
        self.bl_log.info("BLAST database is making")
        os.system('makeblastdb -dbtype nucl -in {0} -out {0}'.format(self.inFile))

        self.bl_log.info("BLAST is running")
        os.system('blastn -query {0} -outfmt 6 -db  {0} -out {1} -window_size 22 -num_threads {2} -evalue {3}'.format(
             self.inFile, self.outFile, self.threads, self.evalue
         ))

        self.filterOut_table()
