import logging
import os
from Bio import SeqIO
from bin.helpers.help_functions import getLog


class TideHunter_run():
    def __init__(self, TH_path, input_fasta, outFasta, THthreads, log_file):
        self.TH_path, self.input_fasta = TH_path, input_fasta
        self.threads = THthreads
        self.outFasta = outFasta
        self.outTab = outFasta + '.tab'
        self.outFasta_all_monomersTH = outFasta + '_all_monomersTH.fasta'
        self.log_th = getLog(log_file, "TideHunter")
        self._run_TH()
        self._tab2fasta()

    def _run_TH(self):
        self.log_th.info("TideHunter has started ....")
        os.system('{0} -f 2 {1} -t {3} > {2}'.format(self.TH_path, self.input_fasta, self.outTab,self.threads))
        self.log_th.info("TideHunter has finished ..... Generation of the consensus sequence is done")
        # os.system('{0} -u {1}  -t {3} > {2}'.format(self.TH_path, self.input_fasta, self.outFasta_all_monomersTH,self.threads))
        # self.log_th.info("TideHunter has finished ..... Generation of the unit sequences of each tandem repeat is done")

    def _tab2fasta(self):
        cnt = 0
        with open(self.outFasta, 'w') as outFasta, open(self.outTab) as outTab:
            for lines in outTab:
                sp = lines.rstrip().split("\t")
                #print(sp)
                seq = ">{0}*{1}*{2}*{3}\n{4}\n".format(sp[0],sp[1], sp[5], sp[6], sp[10])
                outFasta.write(seq)
                cnt += 1
        logging.info("Number of tandem repeats found by TideHunter is: {}".format(cnt))
    
    def getTH2_dic(self):
        """
        return: clusterN:['contig:monomer_seq', ...]
        """
        TH2_fasta = self.outFasta
        dic_to_ret = {} #clustN:[contigN:monomers, ]
        for seq in SeqIO.parse(TH2_fasta, 'fasta'):
            sp = seq.id.split('*')
            if 'artef' not in sp[0]:
                cluster, contig = sp[0].split('_')
                repsNum, monomerLen = sp[-1], sp[-2]
                if cluster not in dic_to_ret:
                    dic_to_ret[cluster] = []
                dic_to_ret[cluster].append(f'{contig}:{str(seq.seq)}')

        return dic_to_ret


