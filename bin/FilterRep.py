import os
import argparse
from Bio import SeqIO
"""
Filtering TR with number of repeats >5 and in second TH file TR monomers for each centroid from louvien clustering table.
"""
class FilteringLouvTab():
    def __init__(self,clustering_outTab,outdir,THall):
        self.clustering_outTab=clustering_outTab
        self.filtering_outTab=outdir+'/louv_clust_filtering.tab'
        self.list_Rep=self.createListRep()
        self.THall_monomers=THall
        self.main(self.list_Rep)
        
        
    def createListRep(self):
        with open(self.clustering_outTab) as LouvTab:
            listFiltRep=[]
            for seqId in LouvTab:
                seq_id=seqId.split('\t')[0]
                numCl=seqId.split('\t')[-1]
                countRep=float(seq_id.split('*')[-1])
                if countRep>5:
                #d4e9d1ee-f7b1-48e7-b8b8-dc7f10aa8435*rep0*368*2.4/0
                    listFiltRep.append('{0}/{1}'.format(seq_id,numCl))
        return listFiltRep    
    
    def main(self,list_Rep):
        with open(self.filtering_outTab,'w') as fastaWr:
            seq_monomer = {}
            for seq in SeqIO.parse(self.THall_monomers, 'fasta'):
                seqID = '*'.join(seq.id.split('_')[0:2])
                if seqID not in seq_monomer:
                    seq_monomer[seqID] = {}
                if seq.id not in seq_monomer[seqID]:
                    seq_monomer[seqID][seq.description] = str(seq.seq)
            seq_list={}
            for repSeq in list_Rep:
                numClust = repSeq.split('/')[-1].rstrip('\n')
                reFormLV = '*'.join(repSeq.split('*')[0:2])
                if reFormLV not in seq_list:
                    seq_list[reFormLV]=repSeq
                    #d4e9d1ee-f7b1-48e7-b8b8-dc7f10aa8435*rep0: d4e9d1ee-f7b1-48e7-b8b8-dc7f10aa8435*rep0*368*2.4/0
            for reSeq in seq_list:
                nClust=seq_list[reSeq].split('/')[-1]
                if reSeq in seq_monomer:
                  
                    for key in seq_monomer[reSeq]:
                        reFSeqMnm = '>{0}/{1}{2}{3}'.format(key, nClust, seq_monomer[reSeq][key],'\n')
                        fastaWr.write(reFSeqMnm)
                        
    
    
