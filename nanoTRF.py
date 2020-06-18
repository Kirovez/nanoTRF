from bin import read_preparation
from bin import run_TideHunter
from bin import Louv_clustering
from bin.helpers.help_functions import checkDir_or_create
from bin.helpers.help_functions import getLog
from bin import run_BLAST

class nanoTRF():
    def __init__(self, reads,  out_directory="test",
                 threads = 4, pathTH = r'/home/ikirov/Tools/TideHunter-v1.4.2/bin/TideHunter'):
        self.outDirectory = checkDir_or_create(out_directory)
        self.reads = reads
        self.log_file = self.outDirectory + '/loging.log'
        LOG = getLog(self.log_file, 'nanoTRF')
        LOG.info("nanoTRF started...")
        self.read_data = ''

        ####### TideHunter parametres ######
        """
        run TH, format tab to fasta ( variable self.outTH_fasta_name). 
        where sequence ids have view as follow: >readName*repN*consLen*copyNum
        """
        self.TH_path = pathTH
        self.threads = threads
        self.outTH_fasta_name = self.outDirectory + "/TH.out.fasta"
        self.TH_data = ''

        ####################################
        ##############BLAST#####################

        self.outFile = self.outDirectory + "/blast.out"
        self.wordsize = 22
        self.evalue = 2
        self.edge_list_after_blast_file = ''
        self.main()



    def main(self):
        self.read_data = read_preparation.PrepareReads(self.reads)
        self.TH_data = run_TideHunter.TideHunter_run(self.TH_path, self.read_data.read_file, self.outTH_fasta_name,
                                                     self.threads, self.log_file)
        self.TH_raw_tab = self.TH_data.outTab

        ##BLAST run###
        blast_module_data = run_BLAST.run_BLAST(self.TH_data.outFasta, self.outFile, self.threads, self.wordsize, self.evalue, self.log_file)
        self.edge_list_after_blast_file = blast_module_data.edge_list_file


        ##Clusetering##
        Louv_clustering.LouvClustering(self.edge_list_after_blast_file, self.outDirectory, self.log_file)




if __name__ == '__main__':
    import sys
    args = sys.argv
    if args[1] == 'test':
        nanoTRF(r'./test_seq/test_seq.fa')
    else:
        # nanoTRF(r'C:\Users\ikirov\Google Диск\Bioinformatics_2018\practice\seq.fastq')
        nanoTRF(sys.argv[1])