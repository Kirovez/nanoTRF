from bin import read_preparation
from bin import run_TideHunter
from bin import Louv_clustering
from bin import FilterRep
from bin import Consensus_Assembly
from bin.helpers.help_functions import checkDir_or_create
from bin.helpers.help_functions import getLog
from bin import run_BLAST
from bin import Delete_dir
import os

class nanoTRF():
    def __init__(self, reads,  out_directory,minAbundancy,consensus_name,threads, pathTH ,log_name,min_overlap, opt_delete):
        self.outDirectory = checkDir_or_create(out_directory)
        self.reads = reads
        self.log_file = self.outDirectory + '/' + log_name
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
        #########################################
        ##############FILTERING##################
        self.minAbundancy=minAbundancy
        ###############CANU#####################
        self.min_overlap = min_overlap
        self.consensus_name = self.outDirectory + '/'+consensus_name
        self.opt_delete = opt_delete
        
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


        ##Clustering##
        louv_module_data=Louv_clustering.LouvClustering(self.edge_list_after_blast_file, self.outDirectory, self.log_file)

        ###Filtering##
        self.clustering_outTab=louv_module_data.clustering_outTab
        Filt_data=FilterRep.FilteringLouvTab(self.clustering_outTab,self.outDirectory,self.reads,self.TH_all_monomers,self.minAbundancy,self.log_file)
        self.tableFilt=Filt_data.filtering_outTab


        ###Canu###
        Consensus_Assembly.ConsAssembly(self.tableFilt,self.outDirectory,self.log_file,self.min_overlap, self.consensus_name)
        self.dir_clust=consensus_out.outdir_clust
        self.dir_canu=consensus_out.outdir_canu
        ###Delete directories###
       
        Delete_dir.Delete_direct(self.dir_clust,self.dir_canu,self.opt_delete,self.log_file)




if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description='A tool to clustering sequences in fasta file and searching  '
                                                 'consensus among the many sequences for each cluster')
    parser.add_argument("-r","--reads", help="Path to FastQ or Fasta file")
    parser.add_argument("-o","--out_directory", help="Path to work directory for output files where will be saved",default='/run_NTRF')
    parser.add_argument("-m","--max_abundancy", help="The proportion of amount lengths all tandem repeats in one cluster to length all the reads",default=0.0001)
    parser.add_argument("-consensus_name",help="File name with consensus sequences, default name is 'consensus.fasta.'",
                        default='consensus.fasta')
    parser.add_argument("-th","--threads", help="Number of threads for running the module Blast", default=4)
    parser.add_argument("-log_file",
                        help="This file list analysis parameters, modules and files, contains messages generated on the various stages of the NanoTRF work. "
                             "It allows tracking events that happens when NanoTRF runs. Default - loging.log",
                        default='loging.log')
    parser.add_argument("-min_Overlap", help="Number of overlapping nucleotides  between repeats in one cluster", default=15)
    parser.add_argument("-path_TH", help="Path to the location of the TideHunter",
                        default=r'/home/ikirov/Tools/TideHunter-v1.4.2/bin/TideHunter')
    parser.add_argument("-opt_delete",help="Remove unncessary large files and directories from working directory", default="d")
    args = parser.parse_args()
    if not os.path.exists(args.reads):
        print("File {} not found!".format(args.reads))
    else:
        print("File {} found...".format(args.reads))
        nanoTRF(args.reads, args.out_directory, minAbundancy=args.max_abundancy, consensus_name=args.consensus_name,
                threads=args.threads, log_name=args.log_file, min_overlap=args.min_Overlap, pathTH=args.path_TH, opt_delete=args.opt_delete)

