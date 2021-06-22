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
from bin import Run_TRF
from bin import Reclustering
from bin import without_TH
import os
import argparse
import os
import shutil

def main():
    args = get_cmdline_args()
    w_TH=args.run_th
    out_trf=args.out_directory
    outDirectory = '{0}/'.format(checkDir_or_create(args.out_directory))
    reads = args.reads
    log_file =outDirectory + args.log_file
    LOG = getLog(log_file, 'nanoTRF')
    LOG.info("nanoTRF started...")

    read_data = ''

    ####### TideHunter parametres ######
    """
    run TH, format tab to fasta ( variable self.outTH_fasta_name). 
    where sequence ids have view as follow: >readName*repN*consLen*copyNum
    """
    TH_path = args.path_TH
    threads = args.threads
    outTH_fasta_name = outDirectory + "TH.out.fasta"
    outFasta_all_monomersTH = ''
    TH_data = ''

    ####################################
    ##############BLAST#####################
    blast_run = args.blast
    makedb = args.makedb
    outFile = outDirectory + "blast.out"
    wordsize = args.wordsize
    evalue = args.evalue
    edge_list_after_blast_file = ''

    ##########################################
    ###############CLUSTERING################
    clustering_outTab = ''
    minCopy= args.min_copy
    ###############CANU#####################
    canu = args.canu
    min_overlap = args.min_Overlap
    consensus_name = outDirectory + args.nano_trf
    tab_name=outDirectory+args.nano_tab

    ###TRF###

    path_TR = args.TRF_run

    ### Reclustering###

    wordsize_f = args.wordsize_f
    perc_abund = args.perc_abund

    ###MAIN###

    ##READ PREPARATION##

    read_data = read_preparation.PrepareReads(reads)
    #########TH##########
    TH_outFasta=''
    TH_raw_tab=''
    TH_all_monomers=''

    if args.run_th:
       run_data=without_TH.without_TH(w_TH[0],outTH_fasta_name,log_file)
       TH_all_monomers=w_TH[1]
       TH_outFasta=run_data.outFasta
    else:
        TH_data = run_TideHunter.TideHunter_run(TH_path, read_data.read_file, outTH_fasta_name,
                                                threads, log_file)
        TH_raw_tab = TH_data.outTab
        TH_all_monomers = TH_data.outFasta_all_monomersTH
        TH_outFasta = TH_data.outFasta

    ##BLAST run###
    blast_module_data = run_BLAST.run_BLAST(blast_run, makedb, TH_outFasta, outFile,
                                            threads, wordsize, evalue, log_file)
    edge_list_after_blast_file = blast_module_data.edge_list_file
    singleton_list = blast_module_data.not_blast

    ##Clustering##
    louv_module_data = Louv_clustering.LouvClustering(edge_list_after_blast_file, outDirectory, log_file)

    ###Filtering##

    clustering_outTab = louv_module_data.clustering_outTab

    Filt_data = FilterRep.FilteringLouvTab(clustering_outTab, singleton_list, outDirectory, reads,
                                           TH_all_monomers, minCopy,log_file)
    tableFilt = Filt_data.filtering_outTab
    abund_tab = Filt_data.clust_abund

    ###Canu###
    consensus_out = Consensus_Assembly.ConsAssembly(canu,tableFilt, singleton_list,
                                                    TH_outFasta, outDirectory, log_file,
                                                    min_overlap,threads)
    dir_clust = consensus_out.outdir_clust
    dir_canu = consensus_out.outdir_canu
    consensus_file=consensus_out.consensus_fasta

    ###TRF###

    TRF_out = Run_TRF.Run_TRF(path_TR, outDirectory, log_file)
    re_blast = TRF_out.dir_trf
    trf_seq = TRF_out.filt_trf

    ###Reclustering###

    reclust_out = Reclustering.Reclustering(consensus_name,tab_name,blast_run, makedb, threads, wordsize_f,trf_seq,
                                            outDirectory, abund_tab, perc_abund, log_file)
    nanoTRF_abund = reclust_out.nanoTRF_abund

    ###Delete directories###
    os.system('rm {0}consensus.fasta*html'.format(out_trf))
    del_log = getLog(log_file, "DELETE")
    for file_t in os.listdir(outDirectory):
        if os.path.isfile(outDirectory+file_t):
            if file_t != outDirectory+'nanoTRF.fasta' or file_t != outDirectory+'TH.out.fasta' or file_t !=outDirectory+'TH.out.fasta.tab' or file_t != outDirectory+'TR_info.tab' or file_t !=outDirectory+ 'loging.log':
                path_t = outDirectory + file_t
                os.remove(path_t)
    for file in os.listdir(args.out_directory):
        if file.startswith('consensus.fasta') and file.endswith('html') and not os.path.isdir(args.out_directory+'/'+file):
            path_t = args.out_directory + file
            os.remove(path_t)
    if  not args.cleanup:
        del_log.info("Removing directories has started...")
        # Delete an entire directory tree - ./clust/, ./canu/ and ./ReBlast/
        shutil.rmtree(dir_canu)
        shutil.rmtree(dir_clust)
        shutil.rmtree(re_blast)
        # Delete an TRF html. reports and unnecessary BLAST files
    
   
    else:
        del_log.info("Directories are not removed")

def get_cmdline_args():

    parser = argparse.ArgumentParser(description='A tool to clustering sequences in fasta file and searching  '
                                                 'consensus among the many sequences for each cluster')
    parser.add_argument("-r", "--reads", help="Path to FastQ or Fasta file")
    parser.add_argument("-pTH", "--path_TH", help="Path to the location of the TideHunter", default='TideHunter')
    parser.add_argument("-T","--run_th", help="  If previously TH was running", nargs=2)
    parser.add_argument("-cu", "--canu", help="Path to the location of the Canu", default='canu')
    parser.add_argument("-trf", "--TRF_run", help="Path to the location of the Tandem Rapeat Finder",default='trf' )
    parser.add_argument("-o", "--out_directory", help="Path to work directory for output files where will be saved")
    parser.add_argument("-b", "--blast", help="Path to blastn executabled", default='blastn')
    parser.add_argument("-mb", "--makedb", help='Path to makeblastdb executable', default='makeblastdb')
    parser.add_argument("-w", "--wordsize", help='Word size for wordfinder algorithm (length of best perfect match)',
                        default=22)
    parser.add_argument("-w_f", "--wordsize_f", help='Word size for Reblusting(length of best perfect match)',
                        default=15)
    parser.add_argument("-ev", "--evalue", help=' Expectation value (E) threshold for saving hits', default=2)
    parser.add_argument("-m", "--min_copy",
                        help="The minimum number of TRs copy in the data",
                        default=100)
    parser.add_argument("-nano", "--nano_trf",
                        help="File name with consensus sequences, default name - nanoTRF.fasta",
                        default='nanoTRF.fasta')
                        
    parser.add_argument("-tab", "--nano_tab",
                        help="Table file with the TRs abundancy ",
                        default='TR_info.tab')                   
                        
                        
    parser.add_argument("-th", "--threads", help="Number of threads for running the module Blast", default=4)
    parser.add_argument("-lg", "---log_file",
                        help="This file list analysis parameters, modules and files, contains messages generated on the various stages of the NanoTRF work. "
                             "It allows tracking events that happens when NanoTRF runs. Default =loging.log",
                        default='loging.log')
    parser.add_argument("-mOVe", "--min_Overlap",
                        help="Number of overlapping nucleotides  between repeats in one cluster", default=10)
    parser.add_argument("-ca", "--perc_abund", help="Minimum value of the TR cluster abundancy", default=0.009)

    parser.add_argument("-c", "--cleanup", default=True, action="store_true",
                        help="Remove unncessary large files and directories from working directory")
    args = parser.parse_args()
    if args.run_th:
        print("nanoTRF running without TideHunter!!!")
    else:
        print("nanoTRF running with TideHunter!!!")
        if not os.path.exists(args.reads):
            print("ERROR!File {} not found!".format(args.reads))
        else:
            print("File {} found...".format(args.reads))
    return args


if __name__ == "__main__":
    main()
