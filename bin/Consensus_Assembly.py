import os
import argparse
from Bio import AlignIO
from Bio import SeqIO

from bin.helpers.help_functions import getLog


class ConsAssembly():
    def __init__(self,canu,filtering_outTab,outdir,log_file,min_overlap,consensus_name):
        self.filtering_outTab=filtering_outTab
        self.min_overlap=min_overlap
        self.consensus_name=consensus_name
        self.outdir=outdir
        self.outdir_clust=outdir+'/clusters/'
        self.outdir_canu=outdir+'/canu/'
        self.canuRun=canu
        self.canu_log = getLog(log_file, "Consensus assembly")
        self.canu_log.info("CONSENSUS ASSEMBLY has started...")
        self.createfile()
        self.dirFile_canu=self.runCanu()        
        self.writeFileCan()

    def createfile(self):
        with open(self.filtering_outTab) as filtFasta:
            if os.path.exists(self.outdir_clust):
                self.canu_log.info("!!ERROR!! Directory specified 'clust' exists.Please check the directory")
            else:
                os.mkdir(self.outdir_clust)
                self.canu_log.info("Directory 'clust' has created")
                self.canu_log.info("Creating file with monomer sequences for each clusters has started...")        
          
                for seq in SeqIO.parse(filtFasta, 'fasta'):
                    name_file = 'clust{}.fasta'.format(seq.id.split('/')[-1])
                    FileFullPath = os.path.join(self.outdir_clust, name_file)
                    if name_file not in os.listdir(self.outdir_clust):
                        t = 'w'
                        self.canu_log.info("File with sequences for {} cluster has created".format(seq.id.split('/')[-1]))
                    else:
                        t = 'a'
                    handle = open(FileFullPath,t)
                    SeqIO.write(seq,handle,'fasta')
                    handle.close()
                self.canu_log.info("Creating all files have finished")

    def runCanu(self):
        if os.path.exists(self.outdir_canu):
            self.canu_log.info("!!ERROR!! Directory specified 'canu' exists.Please check the directory")
        else:
            os.mkdir(self.outdir_canu)
            self.canu_log.info("Directory 'canu' has created")          
          
            for fasta in os.listdir(self.outdir_clust):
                name_Dir=fasta.split('.fasta')[0]
                path_Dir=self.outdir_canu+name_Dir
                os.mkdir(path_Dir)
                run_Canu = '{0} -p {1} -d {2} useGrid=0 -nanopore-raw {3} genomeSize=1000 minReadLength=50 minOverlapLength={4} corMinCoverage=3 stopOnLowCoverage=0'.format(self.canuRun,name_Dir,path_Dir,self.outdir_clust+fasta,self.min_overlap)
                self.canu_log.info("Canu for {} has started".format(name_Dir))
                run = os.system(run_Canu)
          


    def writeFileCan(self):
        self.canu_log.info("Generation consensus file....")
        with open(self.consensus_name,'w') as ConsensusFile:
            countDir = 0
            for dirClust in os.listdir(self.outdir_canu):
                countDir+=1
                path_canudir=self.outdir_canu+dirClust+'/'
                for fasta in os.listdir(path_canudir):
                    fasta_name=dirClust+'.contigs.fasta'
                    if os.path.exists(path_canudir+fasta_name)  :
                        if fasta.endswith('.contigs.fasta'):
                            len_max=0
                            name_tig=''
                            for seq in SeqIO.parse(self.outdir_canu+dirClust+'/'+fasta,'fasta'):
                                sp0 = seq.description.split(' ')
                                len_reads=float(sp0[2].split('=')[-1])
                                if len_reads>len_max:
                                    len_max =len_reads
                                    tig_form='>cons_{2}/{0}_{1}'.format(sp0[1].split('=')[1],sp0[2].split('=')[1],fasta.split('.contig')[0])
                                    seq_p = '{0}\n{1}\n'.format(tig_form,seq.seq)
                                    name_tig=seq_p
                            ConsensusFile.write(name_tig)
                            self.canu_log.info("Consensus for {0} was added in consensus file".format(fasta))
                    else:
                        self.canu_log.info("!!!!Consensus sequences for {0} wasn't assembly!!!!".format(dirClust))
        self.canu_log.info("FINISHED.{} consensus sequences was generated".format(countDir))


 
