"""
makeblastdb -dbtype nucl -in nanoTRF_0.1M.fasta -out nanoTRF_0.1M.fasta

blastn -query merged_TR_rank_all.fasta -outfmt 6 -db  nanoTRF_0.1M.fasta -out merged_TR_rank_all_vs_nanoTRF_0.1M.out -window_size 22 -num_threads 100 -evalue 10
"""
from Bio import SeqIO
from bin.helpers.help_functions import getLog
import os

class run_BLAST():
    def __init__(self,blast_run,makedb,inFile, outFile, threads, wordsize, evalue, log_file, identity_allowed = 90.0, query_sbj_length_differences_allowed = 0.8):
        self.blast_run,self.makedb,self.inFile, self.outFile, self.threads, self.wordsize, self.evalue = blast_run,makedb,inFile, outFile, threads, wordsize, evalue
        self.bl_log = getLog(log_file, "BLAST module")
        self.not_blast=outFile+"_notBlast.list"
        self.edge_list_file = outFile + "edges.list"
        self.query_sbj_length_differences_allowed = query_sbj_length_differences_allowed # ratio between the longest and shortest query and subject in a pair after BLAST.
        self.identity_allowed = identity_allowed
        self.main()

    def filterOut_table(self):
        dict_blast={}
        """
        :return: list of edges [(query, hit),(),...]
        """
        edge_cnt = 0
        singleton_r=0
        count_r=0
        dict_added_edges = {}
        with open(self.outFile) as inFile, open(self.edge_list_file, 'w') as outEdgeList, open(self.not_blast,'w') as outNotBlast:
            for lines in inFile:
                sp = lines.split("\t")
                qseqid, sseqid, evalue, qlen, slen, length, pident = sp
                qlen, slen, pident, length = int(qlen), int(slen), float(pident), int(length)
                len_q_s = [qlen, slen]
                if qseqid not in dict_added_edges:
                    dict_added_edges[qseqid] = []
                if sseqid not in dict_added_edges[qseqid]:
                    dict_added_edges[qseqid].append(sseqid)
                ## FILTERING THE BLAST RESULTS
                    if sp[0] != sp[1] and float(sp[2]) < 0.00001 and (length/qlen) >= self.query_sbj_length_differences_allowed and pident >= self.identity_allowed:
                        outEdgeList.write("{0}\t{1}\t{2}\n".format(qseqid, sseqid, evalue))
                        edge_cnt += 1                   
                        if sp[0] not in dict_blast:
                            dict_blast[sp[0]]=0
            
            for seq in SeqIO.parse(self.inFile,'fasta'):
                if seq.description not in dict_blast:
                    singleton_r+=1
                    count_r+=1
                    outNotBlast.write('{0}\tartef{1}'.format(seq.description,count_r)+'\n')
                    
                    
                
            
                
                    
        print('Number of singletons:{}'.format(singleton_r))
        print("Number of edges", edge_cnt)
        self.bl_log.info("NUmber of edges: {}".format(edge_cnt))


    def main(self):
        self.bl_log.info("BLAST database is making")
        os.system('{0} -dbtype nucl -in {1} -out {1}'.format(self.makedb,self.inFile))

        self.bl_log.info("BLAST is running")
        print('{0} -query {1} -outfmt "6 qseqid sseqid evalue qlen slen length pident" -db  {1} -out {2} -max_hsps 1 -word_size {3} -num_threads {4} -evalue {5}'.format(
             self.blast_run,self.inFile, self.outFile, self.wordsize, self.threads, self.evalue
         ))
        os.system('{0} -query {1} -outfmt "6 qseqid sseqid evalue qlen slen length pident" -db  {1} -out {2} -max_hsps 1 -word_size {3} -num_threads {4} -evalue {5}'.format(
             self.blast_run,self.inFile, self.outFile, self.wordsize, self.threads, self.evalue
         ))

        self.filterOut_table()
