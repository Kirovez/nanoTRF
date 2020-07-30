import os
import argparse
from Bio import SeqIO
from collections import defaultdict
"""
Filtering TR with number of repeats >5 and in second TH file TR monomers for each centroid from louvien clustering table.
"""
class FilteringLouvTab():
    def __init__(self,clustering_outTab,outdir,THall, minAbundancy,log_file):
        self.minAbundancy = minAbundancy
        self.clustering_outTab=clustering_outTab
        self.filtering_outTab=outdir+'/louv_clust_filtering.tab'
        self.filt_log = getLog(log_file, "Filtering")
        self.filt_log.info("Filtering and preparing file with monomer sequences has started...")
        self.list_Rep=self.createListRep()
        self.THall_monomers=THall
        self.main(self.list_Rep)
        
        
    def createListRep(self):
        listFiltRep={} 
        dictRep={}
        cluster_abundancy = defaultdict(int)
        with open(self.clustering_outTab) as LouvTab:
            self.filt_log.info("Filtration on the basis of repeat number for each tandem repeat in genome( >5 rep.)...")
            for seqId in LouvTab: ### seq_id*rep0*len*nrepeats \t cluster 
                 if not seqId.startswith('Sequence'):
                    sp = seqId.rstrip().split('\t')
                    seq_id=sp[0]
                    numCl=sp[-1]
                    countRep=float(seq_id.split('*')[-1])
                    countLen=float(seq_id.split('*')[-2])
                    abund_seq=countRep * countLen
                    len_reads+=abund_seq
                    #numCl:abundancy - 45:556455
                    cluster_abundancy[numCl] += abund_seq                      
                    if countRep>5:
                #d4e9d1ee-f7b1-48e7-b8b8-dc7f10aa8435*rep0*368*2.4/0
                        listFiltRep[seq_id] = numCl
                     else:
                        self.filt_log.info("Repeat number of the{0} less than 5, doesn't come into the further analysis".format(seq_id.rstrip()))                        
        print('Length reads is {}'.format(len_reads))
        self.filt_log.info("Selection of high-copy clusters(with summary length of the tandem repeats in cluster > 100 thousand nucleotides)...")     
        for i in listFiltRep:
            if cluster_abundancy[listFiltRep[i]] / float(len_reads) > 0.0001:
                dictRep = {'*'.join(i.split('*')[0:2]): listFiltRep[i]}
            else:
                self.filt_log.info("{0} has abundancy less than 0,01 % of genome. {1} doesn't come into the further analysis".format(listFiltRep[i], i))          
        return dictRep
   
        def main(self,list_Rep):
            with open(self.filtering_outTab,'w') as fastaWr:
                seq_monomer = {}
                for seq in SeqIO.parse(self.THall_monomers, 'fasta'):
                    seqID = '*'.join(seq.id.split('_')[0:2])             
                    if seqID not in seq_monomer:
                        seq_monomer[seqID] = defaultdict(str)
                    seq_monomer[seqID][seq.description] = str(seq.seq)       
                self.filt_log.info('Selection monomer sequences for tandem repeats after filtering has started......')  
                for reSeq in list_Rep:
                    nClust=list_Rep[reSeq]
                    if reSeq in seq_monomer:
                        for key in seq_monomer[reSeq]:                   
                            reFSeqMnm = '>{0}/{1}\n{2}\n'.format(key, nClust, seq_monomer[reSeq][key])              
                            fastaWr.write(reFSeqMnm)
            self.filt_log.info('Filtering and preparing file with monomer sequences has finished')

