import os
import argparse
from Bio import AlignIO
from Bio import SeqIO

class ConsAssembly():
    def __init__(self,filtering_outTab,outdir):
        self.filtering_outTab=filtering_outTab
        self.outdir=outdir
        self.outdir_clust=outdir+'/clusters/'
        self.outdir_canu=outdir+'/canu/'
        self.canuRun="/home/liza/Tools/canu-2.0/Linux-amd64/bin/canu"
        self.createfile()
        self.dirFile_canu=self.runCanu()
        self.consensus=self.outdir+'/consensus.fasta'
        self.writeFileCan()

    def createfile(self):
        with open(self.filtering_outTab) as filtFasta:
            os.mkdir(self.outdir_clust)
            for seq in SeqIO.parse(filtFasta, 'fasta'):
                name_file = 'clust{}.fasta'.format(seq.id.split('/')[-1])
                FileFullPath = os.path.join(self.outdir_clust, name_file)
                if name_file not in os.listdir(self.outdir_clust):
                    t = 'w'
                else:
                    t = 'a'
                handle = open(FileFullPath,t)
                SeqIO.write(seq,handle,'fasta')
                handle.close()

    def runCanu(self):
        dir_canu=[]
        os.mkdir(self.outdir_canu)
        for fasta in os.listdir(self.outdir_clust):
            name_Dir=fasta.split('.fasta')[0]
            path_Dir=self.outdir_canu+fasta.split('.fasta')[0]
            os.mkdir(path_Dir)
            run_Canu = '{0} -p {1} -d {2} useGrid=0 -nanopore-raw {3} genomeSize=1000 minReadLength=50 minOverlapLength=10 corMinCoverage=3 stopOnLowCoverage=0'.format(self.canuRun,name_Dir,path_Dir,self.outdir_clust+fasta)
            run = os.system(run_Canu)
            dir_canu.append(path_Dir)


   def writeFileCan(self):
        contigs_list=[]
        with open(self.consensus,'w') as ConsensusFile:
            for dirClust in os.listdir(self.outdir_canu):
                for fasta in os.listdir(self.outdir_canu+dirClust+'/'):
                    if fasta.endswith('.contigs.fasta'):
                        len_max=0
                        name_tig=''                    
                        for seq in SeqIO.parse(self.outdir_canu+dirClust+'/'+fasta,'fasta'):                  
                            len_reads= float(seq.description.split(' ')[2].split('=')[-1])
                            if len_reads>len_max:
                                len_max =len_reads                                                                           
                                seq_id = '>{0}/{1}'.format(seq.description, fasta.split('.contig')[0])
                                seq_p = '{0}\n{1}\n'.format(seq_id,seq.seq)
                                name_tig=seq_p
                        ConsensusFile.write(name_tig)
 
