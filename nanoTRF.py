from bin import read_preparation
from bin import run_TideHunter
from bin import Louv_clustering
from bin import FilterRep
from bin import Consensus_Assembly
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
        self.outFasta_all_monomersTH = ''
        self.TH_data = ''

        ####################################
        ##############BLAST#####################

        self.outFile = self.outDirectory + "/blast.out"
        self.wordsize = 22
        self.evalue = 2
        self.edge_list_after_blast_file = ''

        ##########################################
        ###############CLUSTERING################
        self.clustering_outTab = ''

        ###MAIN###
        self.main()



    def main(self):

        ##READ PREPARATION##
        self.read_data = read_preparation.PrepareReads(self.reads)

        #########TH##########
        self.TH_data = run_TideHunter.TideHunter_run(self.TH_path, self.read_data.read_file, self.outTH_fasta_name,
                                                     self.threads, self.log_file)
        self.TH_raw_tab = self.TH_data.outTab
        self.TH_all_monomers=self.TH_data.outFasta_all_monomersTH

        ##BLAST run###
        blast_module_data = run_BLAST.run_BLAST(self.TH_data.outFasta, self.outFile, self.threads, self.wordsize, self.evalue, self.log_file)
        self.edge_list_after_blast_file = blast_module_data.edge_list_file


        ##Clusetering##
        louv_module_data=Louv_clustering.LouvClustering(self.edge_list_after_blast_file, self.outDirectory, self.log_file)

        ###Filtering##
        self.clustering_outTab=louv_module_data.clustering_outTab
        self.Filt_data=FilterRep.FilteringLouvTab(self.clustering_outTab,self.outDirectory,self.TH_all_monomers)
        self.tableFilt=Filt_data.filtering_outTab


        ###Canu###
        Consensus_Assembly.ConsAssembly(self.tableFilt,self.outDirectory)




if __name__ == "__main__":
    import argparse
    import os.path
    parser = argparse.ArgumentParser(description='A tool to clustering sequences in fasta file and searching  '
                                                 'consensus among the many sequences for each cluster')
    parser.add_argument("reads", help="Path to FastQ or Fasta file")
    parser.add_argument("out_directory", help="Path to work directory for output files where will be saved")
    #parser.add_argument("-path_TH",  help="Path to the location of the TideHunter")
    args = parser.parse_args()
    if not os.path.exists(args.reads):
        print("File {} not found!".format(args.reads))
    else:
        print("File {} found...".format(args.reads))
        nanoTRF(args.reads, args.out_directory)
