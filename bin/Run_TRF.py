"""

Run_TRF module - module to locate and display TRs in sequences with using TRF
Input arguments:TRF,consensus_name,outdir,log_file

- TRF- path to the TRF
- consensus_name - path to file with the TRs consensus
- outdir - working directory
- path to the log file

1. def createdir - creating the directory 'ReBlast' for following analysis

2. def TRF - running TRF and selecting the sequences from txt.html report
 Running TRF:
'trf consensus.fasta 2 7 7 80 10 50 500'

Selection of all consensus sequences from the TRF reports which are written into the 'TRF_seq.fasta'

3. def filt_tr
Selection of the longest consensus sequence for each of the TRs and record consensus sequences to the 'seqFilt_trf.fasta'

"""

import os
import linecache
from Bio import SeqIO
from bin.helpers.help_functions import getLog

class Run_TRF():
    def __init__(self,TRF,consensus_name,outdir, log_name):
        self.outdir=outdir
        self.dir_trf=outdir+'/ReBlast/'
        self.run_TRF=TRF        
        self.consensus_name=consensus_name
        self.file_num = self.dir_trf+'/TRF_seq_dr.fasta'
        self.filt_trf = self.dir_trf+ '/seqFilt_trf.fasta'
        self.TRF_log=getLog(self.log_name,'TRF')
        self.TRF_log.info("Module Run_TRF has started the job...")
        self.createdir()
        self.TRF()
        self.filt_tr()
     def createdir(self):  
#creating the directory 'ReBlast' for following analysis      
        if os.path.exists(self.dir_trf):
            self.TRF_log.info("!! Directory specified 'ReBlast' exists !!")
        else:
            os.mkdir(self.dir_trf)
            self.TRF_log.info("Directory 'ReBlast' has created")
            self.TRF_log.info("Creating files with consensus sequences for each monomer has started..")
    def TRF(self): 
        #running TRF and selecting the sequences from txt.html report
        self.TRF_log.info("Creating files with consensus sequences for each monomer has started..")
        #running TRF
        run_trf='{0} {1} 2 7 7 80 10 50 500'.format(self.run_TRF,self.consensus_name)
        self.TRF_log.info('TRF has started')
        os.system(run_trf)
        self.TRF_log.info('TRF has finished')
        self.TRF_log.info('Generation file with TRF consensus pattern has started ...')
        #selection from working direcory all txt.html files - reports for each TRs after running TRF
        with open(self.file_num,'w') as wfile:
            for trf_f in os.listdir(self.outdir):
                if trf_f.endswith('txt.html'):
                    with open(trf_f) as file:
                        seq_id=0
                        for line in file:
                            if line.startswith('Sequence'):
                                seq_id=line.split(' ')[-1]
                            if 'pattern' in line:
                                #cons_clust7/1805_67_358
                                wfile.write(str('>{0}_{1}\n'.format(seq_id.rstrip(),line.split('(')[-1].split(')')[0])))                           
                            if line.startswith('A') or line.startswith('T') or line.startswith('G') or line.startswith('C'):
                                if (line.endswith('A\n') or line.endswith('T\n') or line.endswith('G\n') or line.endswith('C\n')):
                                    wfile.write(str(line.rstrip()+'\n'))
    def filt_tr(self):
        #selection of the longest consensus sequence for each of the TRs
        with open(self.filt_trf,'w') as wfile:
            seq_cons = {}
            list_cons=[]
            for seq in SeqIO.parse(self.file_num,'fasta'):
                
                seq_id='cons_{0}_{1}'.format(seq.id.split('_')[1], seq.id.split('_')[2])
                if seq_id not in seq_cons:
                    seq_cons[seq_id] = []
                    seq_cons[seq_id].append(str(seq.seq))
                else:
                    seq_cons[seq_id].append(str(seq.seq))
            for seq_id in seq_cons:
                count=0
                seq_n=0
                for el_seq in seq_cons[seq_id]:
                    if count==0:
                        seq_n=el_seq
                        count=len(seq_n)
                    else:
                        if len(el_seq)<count:
                            count=len(el_seq)
                            seq_n=el_seq
                wfile.write('>{0}_{2}\n{1}\n'.format(str(seq_id), str(seq_n), count))
                list_cons.append(str(seq_id))
        self.TRF_log.info('Generation file with TRF consensus pattern has finished')
        self.TRF_log.info('Addition sequences not included in the TRF analysis to the output file')
        for seq in SeqIO.parse(self.consensus_name,'fasta'):
            if seq.id not in list_cons:
                with open(self.filt_trf, 'a') as file_all:
                    file_all.write('>{0}\n{1}\n'.format(seq.id,seq.seq))
        self.TRF_log.info('Module Run_TRF has finished the job')
