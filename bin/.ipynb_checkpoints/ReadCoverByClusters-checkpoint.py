from Bio import SeqIO
import numpy as np
import os
import matplotlib.pyplot as plt

class ReadCoverageCalculator():
    def __init__(self,conns_fasta, reads_fasta, outDir, blast = 'blastn', word_size = 16, num_threads = 100, coverage = 0.0, evalue = 0.0000001, pident_cut = 70):
        self.conns_fasta = conns_fasta
        self.reads_fasta = reads_fasta
        self.outDir = outDir
        self.blast_out = f'{outDir}/TRs_vs_reads.tab'
        self.blast = blast
        self.word_size = word_size
        self.num_threads = num_threads
        self.coverage = coverage
        self.pident_cut = pident_cut
        self.evalue = evalue
        self.per_read_coverage, self.bp_masked, self.genAbun_dic = {}, {}, {}
        self.total_read_len = 0
        self.main()


    def _getReadLen(self):
        read_len_zero = {}
        for seq in SeqIO.parse(self.reads_fasta, 'fasta'):
            rlen = len(seq.seq)
            read_len_zero[seq.id] = np.zeros((rlen,), dtype=int)
            self.total_read_len += len(seq.seq)
        return read_len_zero
    
    def getGenomAbund(self):
        genAbun_dic = {} # clusterN: genome abundancy%
        for clusters in self.bp_masked:
            genAbun_dic[clusters] = round((self.bp_masked[clusters] * 100)/self.total_read_len, 3)
        return genAbun_dic
    
    def runBLAST(self):
        os.system('{0} -dbtype nucl -in {1} -out {1}'.format('makeblastdb', self.reads_fasta,self.reads_fasta))
        os.system('{0} -query {1} -outfmt "6 qseqid sseqid sstart send evalue qlen slen length pident" -db  {2} -out {3} -num_alignments 5000 -word_size {4} -num_threads {5}'.format(
                         self.blast, self.conns_fasta, self.reads_fasta, self.blast_out, self.word_size,  self.num_threads
                     ))
        
    
    def parseBlast_Tab(self, blastTab, zero_vec_read_dic):
        with open(blastTab) as tab:
            masked_zero_vectors = {}
            for lines in tab:
                sp = lines.split('\t')
                qseqid, read_id, sstart, send, evalue, qlen, slen, length, pident = sp
                sstart, send, length, pident, qlen, evalue = int(sstart), int(send), int(length), float(pident),int(qlen), float(evalue)
                sstart, send = min(sstart, send), max(sstart, send)

                if length/qlen >= self.coverage and pident >= self.pident_cut and evalue < self.evalue:
                    cluster_id = qseqid.split('_')[0]
                    if cluster_id not in masked_zero_vectors:
                        masked_zero_vectors[cluster_id] = {}
                    if read_id not in masked_zero_vectors[cluster_id]:
                        masked_zero_vectors[cluster_id][read_id] = zero_vec_read_dic[read_id]
                    masked_zero_vectors[cluster_id][read_id][sstart - 1: send] += 1
            return masked_zero_vectors

    def getCove_masked_dics(self, masked_zero_vectors):
        per_read_coverage = {} # cluster: read_coverage%
        bp_masked = {} # cluster: bp masked in all reads
        for clusters in masked_zero_vectors:
            if clusters not in per_read_coverage:
                per_read_coverage[clusters] = []
                bp_masked[clusters] = 0
                
            for read in masked_zero_vectors[clusters]:
                blast_pos_bp = np.count_nonzero(masked_zero_vectors[clusters][read])
                cover = (blast_pos_bp * 100) / len(masked_zero_vectors[clusters][read])
                # if cover > 99:
                #     print(clusters, read)
                per_read_coverage[clusters].append(round(cover))
                bp_masked[clusters] += blast_pos_bp
        return [per_read_coverage, bp_masked]

    def drawHisto(self):
        for clusters in self.per_read_coverage:
            if 'artef' not in clusters:
                plt.style.use('ggplot')
                plt.rcParams["figure.figsize"] = [12, 8]
                font = {'family' : 'sans',
                    'weight' : 'bold',
                    'size'   : 25}

                plt.rc('font', **font)

                plt.title(f"Per_read coverage by {clusters}, %")
                n_bins = 10
                # Plot the histogram
                plt.hist(self.per_read_coverage[clusters], bins = 10)
                # ticklabels = [i for i in range(100, 10)]
                # plt.xticks(ticklabels, ticklabels)

                # Save the histogram
                plt.savefig(f'{self.outDir}/clusters/{clusters}.hist.png')
                plt.clf()
                
    def drawPie(self):
        for clusters in self.per_read_coverage:
            if 'artef' not in clusters:
                plt.style.use('ggplot')
                plt.rcParams["figure.figsize"] = [12, 8]
                font = {'family' : 'sans',
                    'size'   : 40}
                plt.title(f"Pie chat of per_read coverage by {clusters}, %")


                # Pie chart, where the slices will be ordered and plotted counter-clockwise:
                labels = '>90%', '0..90%'
                conv_higher90 = len([i for i in self.per_read_coverage[clusters] if i>=90])
                conv_less90 = len([i for i in self.per_read_coverage[clusters] if i<90])
                sizes = [conv_higher90, conv_less90]

                fig1, ax1 = plt.subplots()
                ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                        shadow=True, startangle=90, textprops=font)
                ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

                # Save the histogram
                plt.savefig(f'{self.outDir}/clusters/{clusters}.pie.png')
                plt.clf()
                
                
    def main(self):
        self.runBLAST()
        zero_vec_read_dic = self._getReadLen()
        masked_zero_vectors = self.parseBlast_Tab(self.blast_out, zero_vec_read_dic)
        self.per_read_coverage, self.bp_masked = self.getCove_masked_dics(masked_zero_vectors)
        self.genAbun_dic = self.getGenomAbund()