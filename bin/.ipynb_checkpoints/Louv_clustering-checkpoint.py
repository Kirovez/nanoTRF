from bin.helpers.help_functions import getLog
import networkx as nx
from community import community_louvain

class LouvClustering():
    def __init__(self, edge_list_file, outdir, log_file):
        self.edge_list_file = edge_list_file
        self.cl_log = getLog(log_file, "Clustering")
        self.cl_log.info("CLUSTERING started ...")
        self.edge_list = self.getEdgeListFromFile()
        self.clustering_outTab = outdir + '/louv_clusering.tab'
        self.G = self.getGraph(self.edge_list)
        self.cl_log.info("Graph partitioning is in process...")
        self.partition = community_louvain.best_partition(self.G)
        self.main()


    def getGraph(self, edge_list):
        self.cl_log.info("Graph generation....")
        print("Graph generation....")
        G = nx.Graph()
        G.add_edges_from(edge_list)
        self.cl_log.info("Graph has been generated from edge list....")
        return G

    def getEdgeListFromFile(self):
        edges = []
        with open(self.edge_list_file) as EdgeFile:
            for lines in EdgeFile:
                sp = lines.rstrip().split("\t")
                edges.append((sp[0], sp[1]))
        return edges
    
    def getInfoClusters(self, clusters_dic):
        cnt_cl_OK = 0
        cnt_zero_cluster = 0
        self.cl_log.info("Community detection is running ...")
        for clusters in clusters_dic:
            if len(clusters_dic[clusters]) > 1:
                cnt_cl_OK +=1
            else:
                cnt_zero_cluster += 1
        self.cl_log.info(f"\n Number of clusters detected: \n with >1 sequences - {cnt_cl_OK} \n with 1 sequence - {cnt_zero_cluster} \n")
        
    def main(self):
        self.cl_log.info("Community detection is running ...")

        clusters_dic = {} #cluster:[seq1, seq2, ...]
        with open(self.clustering_outTab, "w") as out:
            out.write("Sequence\tCluster\n")
            for num, sequence in enumerate(self.partition):
                out.write("{0}\t{1}\n".format(sequence, self.partition[sequence]))
                if self.partition[sequence] not in clusters_dic:
                    clusters_dic[self.partition[sequence]] = []
                clusters_dic[self.partition[sequence]].append(sequence)
        
        self.getInfoClusters(clusters_dic)