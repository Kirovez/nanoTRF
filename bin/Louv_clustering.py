from bin.helpers.help_functions import getLog
import networkx as nx
#import community as community_louvain
from community import community_louvain
class LouvClustering():
    def __init__(self, edge_list_file, outdir, log_file):
        self.edge_list_file = edge_list_file
        self.cl_log = getLog(log_file, "Clustering")
        self.cl_log.info("CLUSTERING started ...")
        self.edge_list = self.getEdgeListFromFile()
        self.clustering_outTab = outdir + '/louv_clusering.tab'
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

    def main(self):
        G = self.getGraph(self.edge_list)
        self.cl_log.info("Community detection is running ...")
        partition = community_louvain.best_partition(G)

        with open(self.clustering_outTab, "w") as out:
            out.write("Sequence\tCluster\n")
            for num, sequence in enumerate(partition):
                out.write("{0}\t{1}\n".format(sequence, partition[sequence]))

