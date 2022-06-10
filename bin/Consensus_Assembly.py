import os
import argparse
from Bio import AlignIO
from Bio import SeqIO
import subprocess
          
import time
from timeit import default_timer as timer
from multiprocessing import Process, Queue, Pool


from bin.helpers.help_functions import getLog


class ConsAssembly():
    def __init__(self,cap3,filtering_outTab,singleton_list,outFasta,outdir,log_file,min_overlap,threads):
        self.threads=threads
        self.filtering_outTab=filtering_outTab
        self.min_overlap=min_overlap
        self.consensus_fasta=outdir+'/consensus.fasta'
        self.outdir=outdir
        self.outFasta=outFasta
        self.outdir_clust=outdir+'/clusters/'
        self.outdir_canu=outdir+'/canu/'
        self.cap3=cap3
        self.dic_assembly_details = {} # clustN:[contig_size min, contig_size max, contig1 sequence]
        self.singleton_list=singleton_list
        self.canu_log = getLog(log_file, "Consensus assembly")
        self.canu_log.info("CONSENSUS ASSEMBLY has started...")
        self.read_fasta_per_cluster_list = self.createfile()
        self.path_to_contigs_all_clusters=self.start_paral_cap3() #self.runCAP3()        
        self.writeFileCAP3(self.path_to_contigs_all_clusters)

    def createfile(self):
        created_files = []
        with open(self.filtering_outTab) as filtFasta:
            if os.path.exists(self.outdir_clust):
                self.canu_log.info("!!ERROR!! Directory specified 'clust' exists.")
            else:
                os.mkdir(self.outdir_clust)
                self.canu_log.info("Directory 'clust' has created")
            self.canu_log.info("Creating file with monomer sequences for each clusters has started...")        
            for seq in SeqIO.parse(filtFasta, 'fasta'):
                name_file = 'clust{}.fasta'.format(seq.id.split('/')[-1])
                FileFullPath = os.path.join(self.outdir_clust, name_file)
                if FileFullPath not in created_files:#os.listdir(self.outdir_clust):
                    t = 'w'
                    self.canu_log.info("File with sequences for {} cluster has created".format(seq.id.split('/')[-1]))
                    created_files.append(FileFullPath)
                else:
                    t = 'a'
                handle = open(FileFullPath,t)
                SeqIO.write(seq,handle,'fasta')
                handle.close()
                
            self.canu_log.info("Creating all files have finished")
        return created_files

    def cap3_run(self,fasta_file, qout):
        f = open(fasta_file + "cap.log" , 'w') 
        run = subprocess.run([self.cap3, fasta_file, '-h', '100', '-n', '-2', '-m', '3', '-p', '80', '-s', '600'],  check=True, stdout=f) #os.system(run_Canu)
        f.close()
        contigs_fasta = fasta_file + '.cap.contigs'
        contigs_ace = fasta_file + '.ace'
        qout.put(contigs_fasta)

    

    def start_paral_cap3(self, threads = 100):
        path_to_contigs_all_clusters = {}
        fasta_files = self.read_fasta_per_cluster_list
        cnt_f = len(fasta_files)
        cnt_0_cap3 = 0 # number of caps file contigs with size 0
        #print('Cluster files:', os.listdir(self.outdir_clust))
        qout = Queue()
        start = timer()

        self.canu_log.info(f'CONSENSUS ASSEMBLY FOR {cnt_f} CLUSTERS IS GOING ....')

        all_processes = []
        for process in range(len(fasta_files)):
            p = Process(target=self.cap3_run, args = (fasta_files[process], qout, ))
            p.start()
            all_processes.append(p)

        for p in all_processes:
            p.join()

        end = timer()
        self.canu_log.info(f'Elapsed time for consensus assembly: {end - start}')
        
        for i in fasta_files:
            contig_generated_by_cap3 = qout.get()
            stats = os.stat(contig_generated_by_cap3)
            #print(contig_generated_by_cap3, stats.st_size, 'bytes')
            clustN = 'clust' + contig_generated_by_cap3.split('clust')[-1].split('.')[0]
            if stats.st_size != 0:
                path_to_contigs_all_clusters[clustN] = contig_generated_by_cap3 #clustN:path to contigs
                #print(path_to_contigs_all_clusters[clustN])
                monomer_size_contig1seq = self.getMonomerSizeContig1seq(path_to_contigs_all_clusters[clustN]) # [min len, max len]
                self.dic_assembly_details[clustN] = monomer_size_contig1seq
            else:
                cnt_0_cap3 += 1
                self.canu_log.warning(f"Cluster {clustN} failed to be assembled by CAP3")
                self.dic_assembly_details[clustN] = ['Cap3 failed','Cap3 failed','Cap3 failed']
        self.canu_log.info(f'Number of clusters with no successful assembly by Cap3:{cnt_0_cap3}')
        return path_to_contigs_all_clusters
    
          
    def getMonomerSizeContig1seq(self, path_to_contigs_all_clusters):
        ms = []
        contig1_seq = ''
        for seq in SeqIO.parse(path_to_contigs_all_clusters, 'fasta'):
            ms.append(len(seq.seq))
            if seq.id == 'Contig1':
                contig1_seq = str(seq.seq)
        return [min(ms), max(ms), contig1_seq]
            
    def writeFileCAP3(self, path_to_contigs_all_clusters):
        self.canu_log.info("Generation consensus file....")
        with open(self.consensus_fasta,'w') as ConsensusFile:
            countDir = 0
            for clust in path_to_contigs_all_clusters:
                countDir+=1
                path_contigs_cap3=path_to_contigs_all_clusters[clust]
                for seq in SeqIO.parse(path_contigs_cap3, 'fasta'):
                    seq.id = clust + "_"  + seq.id
                    SeqIO.write(seq, ConsensusFile, 'fasta')
        self.canu_log.info(f"FINISHED. A single file {self.consensus_fasta} with all contigs was generated") 
        
 