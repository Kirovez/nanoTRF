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
from bin.helpers.help_functions import getRExDB
from bin import run_BLAST
from bin import Run_TRF
from bin import Reclustering
from bin import without_TH
from bin import drawClusters
from bin import ReadCoverByClusters
from bin import diamondAnnotation
from Bio import SeqIO
import os
import argparse
import os
import shutil
from bin import HtmlWriter

def main():
    script_directory = os.path.dirname(__file__)
    print(script_directory)
    args = get_cmdline_args()
    w_TH=args.run_th
    out_trf=args.out_directory
    outDirectory = '{0}/'.format(checkDir_or_create(args.out_directory))
    print(outDirectory)
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
    identity_allowed = args.min_id
    query_sbj_length_differences_allowed = args.query_sbj_length_differences_allowed
    ##########################################
    ###############CLUSTERING################
    clustering_outTab = ''
    minCopy= args.min_copy
    ###############CANU#####################
    cap3 = args.cap3
    min_overlap = args.min_Overlap
    consensus_name = outDirectory + args.nano_trf
    tab_name=outDirectory+args.nano_tab
    min_abundancy_to_draw = args.min_abundancy_to_draw
    mask_blast_word_size = args.mask_blast_word_size  
    mask_blast_query_coverage = args.mask_blast_query_coverage
    mask_blast_identity = args.mask_blast_identity

    
    ### DAIMOND ###
    diamond_path = args.diamond
    if not args.rexdb_fasta:
        rexdb_fasta, rexdb_tab = getRExDB()
    else:
        rexdb_fasta = args.rexdb_fasta
        rexdb_tab = args.rexdb_tab
        

    
    ############################################################################################################################################################
    ####################################################### MAIN ###############################################################################################
    ############################################################################################################################################################
    
    ## 1. READ PREPARATION##

    read_data = read_preparation.PrepareReads(reads)
    
    ## 2. TH ##
    TH_outFasta=''
    TH_raw_tab=''
    TH_all_monomers=''

    if args.run_th:
        run_data=without_TH.without_TH(w_TH,outTH_fasta_name,log_file)
        #TH_all_monomers=w_TH[1]
        TH_outFasta=run_data.outFasta
    else:
        TH_data = run_TideHunter.TideHunter_run(TH_path, read_data.read_file, outTH_fasta_name,
                                                threads, log_file)
        TH_raw_tab = TH_data.outTab
        #TH_all_monomers = TH_data.outFasta_all_monomersTH
        TH_outFasta = TH_data.outFasta

    ## 3. BLAST run ##
    blast_module_data = run_BLAST.run_BLAST(blast_run, makedb, TH_outFasta, outFile,
                                            threads, wordsize, evalue, log_file, identity_allowed = identity_allowed, query_sbj_length_differences_allowed=query_sbj_length_differences_allowed)
    edge_list_after_blast_file = blast_module_data.edge_list_file
    singleton_list = blast_module_data.not_blast

    ## 4. Clustering ##
    louv_module_data = Louv_clustering.LouvClustering(edge_list_after_blast_file, outDirectory, log_file)


    ## 5. Filtering ##
    clustering_outTab = louv_module_data.clustering_outTab

    Filt_data = FilterRep.FilteringLouvTab(clustering_outTab, singleton_list, outDirectory, read_data.read_file,
                                           TH_outFasta , minCopy,log_file) ##### TH_all_monomers
    tableFilt = Filt_data.filtering_outTab
    abund_tab = Filt_data.clust_abund

    ###CAP3###
    consensus_out = Consensus_Assembly.ConsAssembly(cap3,tableFilt, singleton_list,
                                                    TH_outFasta, outDirectory, log_file,
                                                    min_overlap,threads)

    consensus_file=consensus_out.consensus_fasta
    
    
        
    # 6 Cluster drawing object
    draw_log = getLog(log_file, "CLUSTER DRAWING")
    draw_log.info("CLUSTER LAYOUT CALCULATION....ß")
    drawCluster_obj = drawClusters.drawClusters(louv_module_data.G, consensus_out.outdir_clust,louv_module_data.partition)
    
    ## adding singleton TRs
    noBLAST_TR_ids = {}
    with open(singleton_list) as inFile:
        for lines in inFile:
            tr_id, artefactId = lines.rstrip().split('\t')
            noBLAST_TR_ids[tr_id] = artefactId
    with open(consensus_file, 'a') as outFile:
        for seq in SeqIO.parse(TH_outFasta, 'fasta'):
            #print(seq.id)
            if seq.id in noBLAST_TR_ids:
                tr_id = seq.id 
                seq.id, seq.description = noBLAST_TR_ids[seq.id], tr_id
                SeqIO.write(seq, outFile, 'fasta')
                
                
    ###### Run Abundacny estimation by BLAST vs reads ###
    LOG.info('Correction of genome abundancy is going ...')
    ReadCoverByClusters_obj = ReadCoverByClusters.ReadCoverageCalculator(consensus_file, read_data.read_file, outDirectory,
                                                                            word_size = mask_blast_word_size, num_threads = threads, 
                                                                         pident_cut = mask_blast_identity)

    
    LOG.info('Per read coverage estimation and histogram drawign are going ...')
    genome_abund_by_BLAST_vs_reads = ReadCoverByClusters_obj.genAbun_dic
    ReadCoverByClusters_obj.drawHisto()
    ReadCoverByClusters_obj.drawPie()
    
    
    ###### DIAMOND ANNOTATION #####
    diamond_log = getLog(log_file, "DIAMOND")
    diamond_log.info("DIAMOND ANNOTATION START")
    
    ## count number of clusters for annotation ##
    clusters_to_annotate = {}
    for clusters in Filt_data.dic_clust_abund:
        if clusters in Filt_data.dic_clust_abund:
                if 'artef' not in clusters and Filt_data.dic_clust_abund[clusters] > min_abundancy_to_draw:
                    clusters_to_annotate[clusters] = 0
                    
    cnt_cl = len(clusters_to_annotate)
    diamond_log.info(f'{cnt_cl} CLUSTERS WILL BE ANNOTATED BY RExDB and DIAMOND')                


    diamond_obj = diamondAnnotation.DiamondRunAndParse(diamond_path, reads, outDirectory, rexdb_fasta, rexdb_tab)
    n_cl = 0
    for cluster in ReadCoverByClusters_obj.genAbun_dic:
        clust_id = cluster.split('clust')
        if cluster in clusters_to_annotate:
            n_cl += 1
            diamond_log.info(f"Annotation of cluster {n_cl} of {cnt_cl}")
            clust_id = clust_id[1]
            di_hits_per_read_per_classTE, numReadsWithHits, setHits = diamond_obj.runDiamond(clust_id)
            diamond_log.info(f"For cluster {cluster} {numReadsWithHits} reads have similarity to RExDB proteins")

    ###### run TideHuneter on assembled contigs ######
    LOG.info("Tandem repeat finder is running ....")
    outTH_fasta_name_TH2 = outDirectory + "/TH2.search.out.fasta"
    TH2 = run_TideHunter.TideHunter_run(TH_path, consensus_file, outTH_fasta_name_TH2,
                                                threads, log_file)
    trf_dic = TH2.getTH2_dic() #Run_TRF.run_TRF(outDirectory) {clusterN:'contigN:monomer',.....}
    ###################
    
    with open(abund_tab, 'w') as finalTab:
        #print(Filt_data.dic_clust_abund)
        Filt_data.dic_clust_abund = dict(sorted(Filt_data.dic_clust_abund.items(), key=lambda item: item[1], reverse = True))
        #print(dict(sorted(Filt_data.dic_clust_abund.items())))
        finalTab.write(f'Cluster.id\tmin.Contig.Cap3.Length\tmax.Contig.Cap3.Length\tGenome.portion\tcorr.Genome.portion\tContig1.sequence\tSubrepeats.seq\tSubrepeats.len\tAnnotation\n')
        for cluster in Filt_data.dic_clust_abund:
            abundancy = Filt_data.dic_clust_abund[cluster]
            if cluster in consensus_out.dic_assembly_details:
                minContig, maxContig, Contig1_seq = consensus_out.dic_assembly_details[cluster]
                if abundancy > min_abundancy_to_draw:
                    drawCluster_obj.draw_clusters(target_cluster = [int(cluster.split('clust')[1])])
                    draw_log.info(f'Cluster {cluster} has been drawn successfully!')
                    
            else:
                minContig, maxContig, Contig1_seq = 0, 0, ''
                
            if cluster in trf_dic:
                trf_seqs = ';'.join(trf_dic[cluster])
                len_monomers = [len(i.split(':')[1]) for i in trf_dic[cluster]]
                #trf_lens = '\n'.join(str(i) for i in len_monomers)
                minContigTrf,maxContigTrf = min(len_monomers), max(len_monomers)
                
            else:
                trf_seqs = 'nf'
                trf_lens = 'nf'
                minContigTrf,maxContigTrf = 'nf', 'nf'
            corr_abundancy = 0
            if cluster in genome_abund_by_BLAST_vs_reads:
                corr_abundancy = genome_abund_by_BLAST_vs_reads[cluster]
            
            if cluster in diamond_obj.diamond_INFO_per_CLUSTER:
                
                total_read_in_cluster = diamond_obj.reads_in_clusters[cluster]
                di_hits_per_read_per_classTE, numReadsWithHits, setHits = diamond_obj.diamond_INFO_per_CLUSTER[cluster]
                percent_annotated = round((numReadsWithHits*100)/total_read_in_cluster,2)
                hits = ';'.join([i for i in setHits])
                annotation = f'Reads: {numReadsWithHits} ({percent_annotated}%)  : {hits}'
            else:
                annotation = 'nf'
            finalTab.write(f'{cluster}\t{minContig}\t{maxContig}\t{abundancy}\t{corr_abundancy}\t{Contig1_seq}\t{trf_seqs}\t{minContigTrf} - {maxContigTrf}\t{annotation}\n')

    ### write html report ###
    draw_log.info('Html file writing...')
    HtmlWriter.HtmlTabWriter(outDirectory)
    
    del_log = getLog(log_file, "DELETE")
    if  not args.cleanup:
        for file in os.listdir(outDirectory):
            if os.path.isfile(outDirectory+file):
                file=file.rstrip()
                if file!= 'nanoTRF.fasta' and file!= 'TH.out.fasta' and file!='TH.out.fasta.tab' and file!= 'TR_info.tab' and file!='loging.log':
                    path_t = outDirectory + file
                    os.remove(path_t)
        for file in os.listdir(args.out_directory):
            if file.startswith('consensus.fasta') and file.endswith('html') and not os.path.isdir(args.out_directory+'/'+file):
                path_t = args.out_directory + file
                os.remove(path_t)

            del_log.info("Removing directories has started...")
            # Delete an entire directory tree - ./clust/, ./canu/ and ./ReBlast/
            shutil.rmtree(dir_canu)
            shutil.rmtree(dir_clust)
            shutil.rmtree(re_blast)
            # Delete an TRF html. reports and unnecessary BLAST files

   
    else:
        del_log.info("Directories are not removed")
    LOG.info("******* nanoTRF HAS BEEN FINISHED SUCCESSFULLY! *******")
    return LOG

def get_cmdline_args():
    parser = argparse.ArgumentParser(description='A tool to clustering sequences in fasta file and searching  '
                                                 'consensus among the many sequences for each cluster')
    parser.add_argument("-r", "--reads", help="Path to FastQ or Fasta file")
    parser.add_argument("-pTH", "--path_TH", help="Path to the location of the TideHunter", default='TideHunter')
    parser.add_argument("-T","--run_th", help="If you do not want to run TideHunter again and you have table file (-f 2 option in Tide Hunter), type the path to this file here")
    parser.add_argument("-cap", "--cap3", help="Path to the location of the Cap3", default='cap3')
    parser.add_argument("-diamond", "--diamond", help="Path to the location of DIAMOND",default='diamond' )
    
    parser.add_argument("-o", "--out_directory", help="Path to work directory for output files where will be saved")
    parser.add_argument("-b", "--blast", help="Path to blastn executabled", default='blastn')
    parser.add_argument("-mb", "--makedb", help='Path to makeblastdb executable', default='makeblastdb')
    parser.add_argument("-w", "--wordsize", help='Word size for wordfinder algorithm (length of best perfect match)',
                        default=24)
    parser.add_argument("-w_f", "--wordsize_f", help='Word size for Reblusting(length of best perfect match)',
                        default=15)
    parser.add_argument("-ev", "--evalue", help=' Expectation value (E) threshold for saving hits', default=2)
    parser.add_argument("-mid", "--min_id", help=' minimum identity between monomers to be selected for clustering', default=80.0, type = float)
    parser.add_argument("-bld", "--query_sbj_length_differences_allowed", help='maximum differences in length between query and subject', default=0.4, type = float )
    parser.add_argument("-mad", "--min_abundancy_to_draw", help = "Minimum genome abundancy for cluster of repeats to be drawn", default= 0.01, type = float)
    
    parser.add_argument("-m", "--min_copy",
                        help="The minimum number of TRs copy in the data",
                        default=100)
    parser.add_argument("-nano", "--nano_trf",
                        help="File name with consensus sequences, default name - nanoTRF.fasta",
                        default='nanoTRF.fasta')
                        
    parser.add_argument("-tab", "--nano_tab",
                        help="Table file with the TRs abundancy ",
                        default='TR_info.tab')  
    
    parser.add_argument("-rexdb_fasta", "--rexdb_fasta",
                        help="Fasta file with the RExDB protein sequences")
                        
    parser.add_argument("-rexdb_tab", "--rexdb_tab",
                        help="Table file with the RExDB annotation")
        
    parser.add_argument("-th", "--threads", help="Number of threads for running the module Blast and TideHunter", default=4)
    parser.add_argument("-lg", "---log_file",
                        help="This file list analysis parameters, modules and files, contains messages generated on the various stages of the NanoTRF work. "
                             "It allows tracking events that happens when NanoTRF runs. Default =loging.log",
                        default='loging.log')
    parser.add_argument("-mOVe", "--min_Overlap",
                        help="Number of overlapping nucleotides  between repeats in one cluster", default=10)
    parser.add_argument("-ca", "--perc_abund", help="Minimum value of the TR cluster abundancy", default=0.009)

    parser.add_argument("-c", "--cleanup", default=True, action="store_false",
                        help="Remove unncessary large files and directories from working directory")
    
    parser.add_argument("-maskws", "--mask_blast_word_size", help='word size of blastn masking of raw reads by cluster contig sequences', default=28, type =int)
    parser.add_argument("-maskcov", "--mask_blast_query_coverage", help='query (contig sequence) coverage in blastn masking of raw reads by cluster contig sequences', default=0.5, type = float)
    parser.add_argument("-maskiden", "--mask_blast_identity", help='minimum identity between query (contig sequence)  and raw reads in blastn masking of raw reads by cluster contig sequences', default=70, type =int)
    
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
    import time
    from timeit import default_timer as timer
    start = timer()    
    LOG = main()
    end = timer()
    print(f'elapsed time: {end - start}')
    LOG.info(f"Elapsed time: {end - start}")
