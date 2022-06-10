import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
from community import community_louvain

class drawClusters():
    def __init__(self, G, outDir, partition):
        self.G = G
        self.pos = nx.spring_layout(self.G)
        self.outDir = outDir
        self.node_size = []
        self.node_color = []
        self.partition = partition
    
    def _getColor(self, repeats):
        dic_colors = {2:'black', 5:'blue', 10: 'green', 20:'orange', 40: 'red'}
        for size in dic_colors:
            if repeats < size:
                return dic_colors[size]
        return 'purple'
            
    def _getTargetNodes(self, target_cluster = []):
        nodes_dic = {}
        for node in self.partition:  #node is a sequence id
            if self.partition[node] in target_cluster:
                nodes_dic[node] = 0
                
        edges_list = []
        self.node_color = [self._getColor(float(i.split('*')[-1])) for i in nodes_dic]
        self.node_size = [(float(i.split('*')[-1])) for i in nodes_dic]
        #print(self.node_color)
        for edges in self.G.edges():
            if edges[0] in nodes_dic or edges[1] in nodes_dic:
                edges_list.append(edges)
        return [[i for i in nodes_dic], edges_list]

    def draw_clusters(self, target_cluster = [53]):
        node_list, edge_list = self._getTargetNodes(target_cluster = target_cluster)
        if node_list:
            font = {'family' : 'sans',
                    'weight' : 'bold',
                    'size'   : 25}

            plt.rc('font', **font)
            print(f'Number of nodes selected is {len(node_list)}')
           #cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
            plt.figure(figsize=(12,8))
            plt.title(f"Cluster clust{'_'.join([str(i) for i in target_cluster])}", fontsize=20)
            # draw the graph
            cmap = plt.cm.get_cmap('rainbow')
            layout = self.pos
            vmin = min(self.node_size)
            vmax = max(self.node_size)

            nx.draw_networkx_edges(self.G, self.pos, alpha=0.5, edgelist = edge_list)
            nod = nx.draw_networkx_nodes(self.G, self.pos, nodelist=node_list,  node_size=100, 
                                   node_color = self.node_size, cmap = cmap, vmin=vmin, vmax=vmax)
            nod.set_edgecolor('darkgrey')
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
            sm.set_array([])
            cbar = plt.colorbar(sm)
            cbar.ax.set_ylabel('Repeats in a read',labelpad=15,rotation=270)
            plt.savefig(f"{self.outDir}/clust{'_'.join([str(i) for i in target_cluster])}.png") 
            #plt.show()
