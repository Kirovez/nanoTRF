"""
NanoTRF is a pipeline for de novo identification and sequence assembly of high-copy tandem repeats in raw data ONT plant DNA data,
to calculate TRs copy number variation and assembly consensus sequences for further analyses (e.g. FISH). 

NanoTRF-based TRs search and reconstruction of the consensus sequences has 9 main steps: 
    1) preparation read for following analysis - converting FASTQ format into FASTA (if input file in FASTQ format) - module read_preparation
    2) tandem repeat detection - module run_TideHunter,
    3) searching for the similarity in previously identifying on the 1st step tandem repeats  - module run_BLAST,
    4) clustering repeat sequences - module Louv_clustering,
    5) selection of the high-copy TRs  - module FiltRep,
    6) producing consensus sequences (5) - module Consensus_Assembly, 
    7) searching for the location and display TRs in sequences  - module Run_TRF,
    8) re-clustering obtained TRs pattern - module Reclustering,
    9) removing unnecessary large files and directories from working directory - module Delete_dir

"""
from bin import read_preparation
from bin import run_TideHunter
from bin import Louv_clustering
from bin import FilterRep
from bin import Consensus_Assembly
from bin.helpers.help_functions import checkDir_or_create
from bin.helpers.help_functions import getLog
from bin import run_BLAST
from bin import Delete_dir
from bin import Run_TRF
from bin import Reclustering
import os

class nanoTRF():
    def __init__(self, reads, path_TH, canu,out_directory,blast, makedb,wordsize,wordsize_f,evalue,minAbundancy,consensus_name,threads,log_name,min_overlap,path_TR opt_delete):
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
        self.TH_path = path_TH
        self.threads = threads
        self.outTH_fasta_name = self.outDirectory + "/TH.out.fasta"
        self.outFasta_all_monomersTH = ''
        self.TH_data = ''

        ####################################
        ##############BLAST#####################
        self.blast_run=blast
        self.makedb=makedb
        self.outFile = self.outDirectory + "/blast.out"
        self.wordsize = wordsize
        self.evalue = evalue
        self.edge_list_after_blast_file = ''

        ##########################################
        ###############CLUSTERING################
        self.clustering_outTab = ''
        self.minAbundancy=minAbundancy
        ###############CANU#####################
        self.canu=canu
        self.min_overlap = min_overlap
        self.consensus_name = self.outDirectory + '/'+consensus_name
        self.opt_delete = opt_delete
        
        ###TRF###
        
        self.path_TR=path_TR
        
        ### Reclustering###
        
        self.wordsize_f=wordsize_f
        
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
        blast_module_data = run_BLAST.run_BLAST(self.blast_run,self.makedb,self.TH_data.outFasta, self.outFile, self.threads, self.wordsize, self.evalue, self.log_file)
        self.edge_list_after_blast_file = blast_module_data.edge_list_file


        ##Clustering##
        louv_module_data=Louv_clustering.LouvClustering(self.edge_list_after_blast_file, self.outDirectory, self.log_file)

        ###Filtering##
        
        self.clustering_outTab=louv_module_data.clustering_outTab
        
        Filt_data=FilterRep.FilteringLouvTab(self.clustering_outTab,self.outDirectory,self.reads,self.TH_all_monomers,self.minAbundancy,self.log_file)
        self.tableFilt=Filt_data.filtering_outTab
        self.abund_tab=Filt_data.clust_abund


        ###Canu###
        consensus_out=Consensus_Assembly.ConsAssembly(self.canu,self.tableFilt,self.outDirectory,self.log_file,self.min_overlap, self.consensus_name)
        self.dir_clust=consensus_out.outdir_clust
        self.dir_canu=consensus_out.outdir_canu
        
        
        ###TRF###
        
        TRF_out=Run_TRF.Run_TRF(self.path_TRF,self.consensus_name,self.log_file)
        self.re_blast=TRF_out.dir_trf
        self.trf_seq=TRF_out.file_trf
        
        
        ###Reclustering###
        
        reclust_out=Reclustering.Reclustering(self.blast_run,self.makedb,self.threads,self.word_trf,self.trf_seq,self.outDirectory,self.abund_tab,self.log_file)
        
        
        ###Delete directories###
       
        Delete_dir.Delete_direct(self.outDirectory,self.dir_clust,self.dir_canu,self.re_blast, self.opt_delete,self.log_file)



if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description='A tool to clustering sequences in fasta file and searching  '
                                                 'consensus among the many sequences for each cluster')
                                      
    parser.add_argument("-r","--reads", help="Path to FastQ or Fasta file")
    parser.add_argument("-pTH", "--path_TH", help="Path to the location of the TideHunter")
    parser.add_argument("-cu","--canu", help="Path to the location of the Canu") 
     parser.add_argument("-trf","--TRF", help="Path to the location of the Tandem Rapeat Finder") 
    parser.add_argument("-out","--out_directory", help="Path to work directory for output files where will be saved")
    parser.add_argument("-bn","--blast", help="Path to blastn executabled",default='blastn')
    parser.add_argument("-mb","--makedb", help='Path to makeblastdb executable', default='makeblastdb')
    parser.add_argument("-w","--wordsize", help='Word size for wordfinder algorithm (length of best perfect match)', default=22)
    parser.add_argument("-w_f","--wordsize_f", help='Word size for Reblusting(length of best perfect match)', default=15)
    parser.add_argument("-ev","--evalue", help=' Expectation value (E) threshold for saving hits', default=2)   
    parser.add_argument("-m","--max_abundancy", help="The proportion of amount lengths all tandem repeats in one cluster to length all the reads",default=0.0001)
    parser.add_argument("-cons","--consensus_name",help="File name with consensus sequences, default name - consensus.fasta",
                        default='consensus.fasta')
    parser.add_argument("-th","--threads", help="Number of threads for running the module Blast", default=4)
    parser.add_argument("-lg","---log_file",
                        help="This file list analysis parameters, modules and files, contains messages generated on the various stages of the NanoTRF work. "
                             "It allows tracking events that happens when NanoTRF runs. Default =loging.log",
                        default='loging.log')
    parser.add_argument("-mOVe","--min_Overlap", help="Number of overlapping nucleotides  between repeats in one cluster", default=15)
  
    parser.add_argument("-del","--opt_delete",help="Remove unncessary large files and directories from working directory", default="d")
    args = parser.parse_args()
    if not os.path.exists(args.reads):
        print("File {} not found!".format(args.reads))
    else:
        print("File {} found...".format(args.reads))
        nanoTRF(reads=args.reads, path_TH=args.path_TH,canu=args.canu, out_directory=args.out_directory,blast=args.blast, makedb=args.makedb,wordsize=args.wordsize,wordsize_f=args.wordsize_f,evalue=args.evalue,minAbundancy=args.max_abundancy, consensus_name=args.consensus_name,
                threads=args.threads, log_name=args.log_file, min_overlap=args.min_Overlap,path_TR=args.TRF,opt_delete=args.opt_delete)
