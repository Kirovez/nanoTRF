"""

'Reclustering module' - clustering searched TRF tandem repeats to group together similar TRs
Input arguments:blast_run, makedb,threads, word_size, trf_file,outdir,abund_f,log_file

- blast_run - path to the BLAST
- makedb - path to the makeblastsdb
- threads - number of threads
- word-size - word length for BLAST analysis
- trf_file - file with TRs after TandemRepeatFinder
- outdir - working diretory
- abunf_f - table with cluster abundancy
- log_file - path to the log file

1. def BLAST -  all-versus-all comparisons of TRs (with using BLASTn)

    makeblastdb -dbtype nucl -in nanoTRF.fasta -out nanoTRF.fasta

    blastn -query nanoTRF.fasta  '-db  nanoTRF.fasta  -out nanoTRF.fasta.out -outfmt'6 qseqid sseqid pident qcovs -window_size 15 -num_threads 100

2. def Blast_parsing -  parsing BLAST table.

    Selection only BLAST pairs that have coverage >=80% and identity>=80%. After filtering BLAST pair are written in the form 'clust_1, clust_2', where the first element is first cluster node,
    second element - next cluster node and clust_1 - clust_2 is cluster edge.
    The BLAST filtering pairs in the list are snown below:
    [clust_1,clust_2
    clust_2,clust_3
    clust_4,clust_6
    ...             ]

3.  def createGraph - graph creating.

 Create an empty graph with no nodes and no edges and add nodes and edges from list with BLAST filtering pairs:

Graph_TR.add_node(edge_f)
Graph_TR.add_node(edge_s)
Graph_TR.add_edge(edge_f, edge_s)

The structure of G is analyzed using nx.connected_components(G) function.
For created graphs is searched for the centroid is the TRs which is linked to the largest number of nodes. All centroid are added to centroid_list.
Also for each graph searched for all nodes related to the centroid.

Centroid and nodes related to the graph center are written into the file 'seq_clust.clst' are shown below:
'cluster_6|cluster_5 '    ,where the first element is the cluster centroid, the second element is the node

4.  def filt_clust - creating fasta file 'TR_nanotrf.fasta'  with the filtering TRS and file 'abund_nanotrf.tab' with cluster abundancy

    Adding centroid TRs and TRs that were not included into Graph-analysis for lack similarity or low similarity with other TRs to created file 'TR_nanotrf.fasta'
    The recalculation of cluster abundancy for TRs cluster  which are similar to one another more then >=80% and adding abundancy values for these clusters.

"""
import os
import networkx as nx
from bin.helpers.help_functions import getLog
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
class Reclustering():
    def __init__(self, blast_run, makedb,threads, word_size, trf_file,outdir,abund_f,log_file):
        self.blast_run,self.makedb,self.outdir,self.threads, self.word_size = blast_run,makedb,outdir, threads, word_size
        self.trf_file=trf_file
        self.outdir_reblast=self.outdir +'/ReBlast/'
        self.out_blast=self.outdir_reblast+'/blast_sec.tab'
        self.out_clust = self.outdir_reblast + '/seq_clust.clst'
        self.nanoTRF=self.outdir+'/TR_nanotrf.fasta'
        self.Reclust_log=getLog(log_file,'Reclustering')
        self.abund_f=abund_f
        self.nanoTRF_abund = self.outdir+'/abund_nanotrf.tab'
        self.BLAST()
        self.list_BLAST =self.Blast_parsing()
        self.fasta_clust=self.createGraph(self.list_BLAST)
        self.filt_clust(self.fasta_clust)

    def BLAST(self):
        #all-versus-all comparisons of TRs
        self.Reclust_log.info('BLAST database  is making')      
        os.system('{0} -in {1} -out seqdb -dbtype nucl'.format(self.makedb, self.trf_file))
        self.Reclust_log.info('Blastn with file {0} is running'.format(self.trf_file))
        os.system("{0} -query {1} -db seqdb -out {2} -outfmt '6 qseqid sseqid pident qcovs' -num_threads {3} -word_size {4}".format(self.blast_run,self.trf_file, self.out_blast, self.threads,self.word_size))
        #filtering  and removing similar strings and sequences with identity <80 and coverage <80
    def Blast_parsing(self):
        #parsing BLAST table
        l_BLAST = []
        self.Reclust_log.info('Table filtration after Blastn is running'+'\n')
        with open(self.out_blast) as pars_tab:
            count=0
            for stB in pars_tab:
                id_q=stB.split('\t')[0]
                id_s=stB.split('\t')[1]
                v_ident=stB.split('\t')[2]
                q_cov=stB.split('\t')[3]
                if id_q!=id_s and float(v_ident)>=80 and float(q_cov)>=80:
                    id_q_s='{0},{1}'.format(id_q,id_s)
                    id_s_q = '{1},{0}'.format(id_q, id_s)
                    if (id_q_s and id_s_q) not in l_BLAST:
                        count+=1
                        l_BLAST.append(id_q_s)
            self.Reclust_log.info('{0} pairs have been added for clustering analysis\n'.format(count))
            self.Reclust_log.info('Table filtration after Blastn is done')
        return l_BLAST

    def createGraph(self,list_BLAST):
        #graph creating
        centroid_name = []
        self.Reclust_log.info('Edge addition is running\n')
        with open(self.out_clust,'w') as clust_t:
            count_clust = 0
            Graph_TR = nx.Graph()
            for edge in list_BLAST:
                edge_f = edge.split(',')[0]
                edge_s = edge.split(',')[-1]
                Graph_TR.add_node(edge_f)
                Graph_TR.add_node(edge_s)
                Graph_TR.add_edge(edge_f, edge_s)
            self.Reclust_log.info('Appending the nodes and edges is done.Generate connected components\n')
            self.Reclust_log.info('Creating graph is running'+'\n')
            connected_components = nx.connected_components(Graph_TR)
            self.Reclust_log.info('Number of clusters: {0}'.format(nx.number_connected_components(Graph_TR)))
            for num, clusters in enumerate(connected_components):
                count_clust+=1
                #Appending the number nodes in cluster
                node_list = []
                degree_centroid=0
                name_centroid=0
                for element in clusters:
                    #number of element degree
                    degree_elem=Graph_TR.degree[element]
                    node_list.append(element)
                    # degree_list.append(degree_elem)
                    if degree_elem > degree_centroid:
                        degree_centroid=degree_elem
                        name_centroid=element
                self.Reclust_log.info('Centroid is {0} for {1} cluster'. format(name_centroid,count_clust)+'\n')
                for num in node_list:
                    if num!=name_centroid:
                        ref_num='cluster_{}'.format(num.split('/')[0].split('clust')[-1])
                        ref_cent= 'cluster_{}'.format(name_centroid.split('/')[0].split('clust')[-1])
                        centroid_name.append(num)
                        clust_t.write('{0}|{1}\n'.format(ref_cent,ref_num))
            self.Reclust_log.info('Searching of centroid in all clusters is done. Creating graph is running\n')
            self.Reclust_log.info('Сentroid search in clusters is done')
            self.Reclust_log.info('Creating graph is done')
        return centroid_name

    def filt_clust(self,fasta_clust):
        elem_tab=[]
        abund_dict={}
        abund_dict2={}       
        with open(self.nanoTRF,'w') as tr_nano:
            for seq in SeqIO.parse(self.trf_file,'fasta'):
                if seq.id not in fasta_clust:
                    SeqIO.write(seq,tr_nano,'fasta')
        #cons_clust48/1949_32_167|cons_clust15/1215_56_167
        with open(self.abund_f) as abund_f, open(self.out_clust) as clust_t:
            for st in abund_f:
                sp=st.rstrip().split('\t')
                if sp[0] not in abund_dict:
                    abund_dict[sp[0]]=sp[1]
            for seq in clust_t:
                sp=seq.split('|')
                sp1 = seq.split('|')[-1].rstrip()
                if sp[0] not in abund_dict2:
                    abund_dict2[sp[0]]=[]
                    abund_dict2[sp[0]].append(abund_dict[sp[0]])
                    abund_dict2[sp[0]].append(abund_dict[sp1])
                    elem_tab.append(sp[1].rstrip())
                else:
                    abund_dict2[sp[0]].append(abund_dict[sp1])
                    elem_tab.append(sp[1].rstrip())
            for seq in abund_dict:
                if seq not in abund_dict2 and seq not in elem_tab:
                    abund_dict2[seq]=abund_dict[seq]
        with open(self.nanoTRF_abund,'w') as nano_abund:
            nano_abund.write('Cluster number\tAbundancy,%\n')
            for seq in abund_dict2:
                count_ab=0
                if type(abund_dict2[seq]) is list:
                    for el in abund_dict2[seq]:
                        count_ab+=float(el)
                else:
                    count_ab+=float(abund_dict2[seq])
                nano_abund.write('{0}\t{1:.5f} \n'.format(seq,count_ab*100))
