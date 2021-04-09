from bin.helpers.help_functions import getLog
import os
import argparse
from Bio import SeqIO
from collections import defaultdict
"""
Filtering TR with number of repeats >5 and in second TH file TR monomers for each centroid from louvien clustering table.
"""
############singletonR
class FilteringLouvTab():
    def __init__(self,clustering_outTab,singleton_list,outdir,reads,THall,minAbundancy,log_file):
        self.minAbundancy = minAbundancy
        self.reads=reads
        self.singletonR=singleton_list
        self.clustering_outTab=clustering_outTab
        self.filtering_outTab=outdir+'/louv_clust_filtering.tab'
        self.clust_abund=outdir+'/clust_abund.tab'
        self.filt_log = getLog(log_file, "Filtering")
        self.filt_log.info("Filtering and preparing file with monomer sequences has started...")
        self.list_Rep=self.createListRep()
        self.THall_monomers=THall
        self.main(self.list_Rep)
        
        
        
    def createListRep(self):
        listFiltRep={}
        dictRep={}
        cluster_abundancy = defaultdict(int)
        copyNum=defaultdict(int)
        len_reads=0
        abun_cl=[]
        list_clust=[]
        for seq in SeqIO.parse(self.reads,'fasta'):
            len_reads+=len(seq.seq)
        with open(self.clustering_outTab) as LouvTab, open(self.singletonR) as singl_R:            
            for seqId in LouvTab: ### seq_id*rep0*len*nrepeats \t cluster 
                 if not seqId.startswith('Sequence'):                   
                    sp = seqId.rstrip().split('\t')
                    seq_id=sp[0]
                    numCl=sp[-1]
                    countRep=float(seq_id.split('*')[-1])
                    countLen=float(seq_id.split('*')[-2])
                    abund_seq=countRep * countLen                  
                    #numCl:abundancy - 45:556455
                    cluster_abundancy[numCl] += abund_seq   
                    copyNum[numCl]+=countRep   
             ##################################################new, version3###################################################
            for seq in singl_R:                
                sp = seq.rstrip().split('\t')
                countRep = float(sp[0].split('*')[-1])
                # append singleton reads and TR copy number in reads
                copyNum[sp[-1]] = countRep
                countLen = float(sp[0].split('*')[-2])
                abund_seq = countRep * countLen
                cluster_abundancy[sp[-1]] += abund_seq      
            for clust in copyNum:
                if copyNum[clust]>=100:
                    list_clust.append(clust)
            print(list_clust)
            
        with open(self.clustering_outTab) as LouvTab, open(self.singletonR) as singl_R:
            for seq in LouvTab:
                if not seq.startswith('Sequence'):
                    sp=seq.rstrip().split('\t')
                    if sp[-1] in list_clust:                        
                        listFiltRep[sp[0]] = sp[-1]  
            for seq in singl_R:
                sp=seq.rstrip().split('\t')
                if sp[-1] in list_clust:
                    listFiltRep[sp[0]] = sp[-1]
                                
        self.filt_log.info('Length reads is {}'.format(len_reads))      
        with open(self.clust_abund,'w') as clust_f:               
            for i in listFiltRep:
                if '*'.join(i.split('*')[0:2]) not in dictRep:
                    dictRep['*'.join(i.split('*')[0:2])]=listFiltRep[i]
                    name_cl='cluster_{0}'.format(listFiltRep[i])
                    if  name_cl not in abun_cl:
                        clust_f.write('cluster_{0}\t{1}\t{2}\n'.format(listFiltRep[i],cluster_abundancy[listFiltRep[i]]/ float(len_reads), cluster_abundancy[listFiltRep[i]]))
                        abun_cl.append(name_cl)
      #dictRep = {'*'.join(i.split('*')[0:2]):listFiltRep[i] for i in listFiltRep if cluster_abundancy[listFiltRep[i]] / float(len_reads) > float(self.minAbundancy)}
        return dictRep
   
    
    def main(self,list_Rep):
        with open(self.filtering_outTab,'w') as fastaWr:
            seq_monomer = {}
            for seq in SeqIO.parse(self.THall_monomers, 'fasta'):
                seqID = '*'.join(seq.id.split('_')[0:2])             
                if seqID not in seq_monomer:
                    seq_monomer[seqID] = {}
                if seq.id not in seq_monomer[seqID]:
                    seq_monomer[seqID][seq.description] = str(seq.seq)       
            self.filt_log.info('Selection monomer sequences for tandem repeats after filtering has started......')                    
            for reSeq in list_Rep:
                nClust=list_Rep[reSeq]
                if 'artef' not in nClust:
                
                    if reSeq in seq_monomer:
                        for key in seq_monomer[reSeq]:                   
                            reFSeqMnm = '>{0}/{1}\n{2}\n'.format(key, nClust, seq_monomer[reSeq][key])              
                            fastaWr.write(reFSeqMnm) 
                        
        self.filt_log.info('Filtering and preparing file with monomer sequences has finished')
