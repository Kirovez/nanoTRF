from collections import defaultdict
from Bio import SeqIO
import os
def _getSeq_per_cluster(clustering_outTab):
    seq_per_cluster = {} #seq_id = cluster
    with open(clustering_outTab) as louv:
        for lines in louv:
            seq_id, cluster_id = lines.rstrip().split("\t")
            seq_per_cluster[seq_id] = cluster_id
    print("Number of clusters", len(seq_per_cluster))
    return seq_per_cluster

def _closeFiles(d_files):
    for f in d_files:
        d_files[f].close()

def collectSeqsFromClusters(outTH_fasta_name, clustering_outTab, outDir, rep_min_singlets = 10000):
    """
    creates multiple fasta file containing monomer sequences for each cluster
    :param louv_tab: path to the tables after Louv_clustering
    :return: path to the folder with fastas
    """
    cluster_per_seq = _getSeq_per_cluster(clustering_outTab)
    max_cluster_num = len(cluster_per_seq) # to add cluster id for singlet TRs
    files = {} #cluster_id:file.fasta

    #iterate through all TR sequences and add to the corresposnding file. If TR not in clustering tabla (=singlet), cluster id will be generated
    for seq in SeqIO.parse(outTH_fasta_name, 'fasta'):
        low_copy = False
        #singlet TRs (not in louv clustering file)
        if seq.id in cluster_per_seq:
            cluster_id = cluster_per_seq[seq.id]
        else:
            sp = seq.id.split("*")
            abundancy = float(sp[-2]) * float(sp[-1])
            if abundancy > rep_min_singlets:
                cluster_id = max_cluster_num + 1
                max_cluster_num += 1
            else:
                low_copy = True

        if not low_copy:
            fasta_file_name = "CL{}".format(cluster_id)

            if fasta_file_name not in files:
                files[fasta_file_name] = open(os.path.join(outDir, "", fasta_file_name), 'w')

            SeqIO.write(seq, files[fasta_file_name], 'fasta')

    _closeFiles(files)
import sys
collectSeqsFromClusters(sys.argv[1], sys.argv[2], sys.argv[3])
