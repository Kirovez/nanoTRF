import os
from Bio import SeqIO
import matplotlib.pyplot as plt

class DiamondRunAndParse():
    def __init__(self, diamond_path, raw_fasta_reads, outDir, rexdb_fasta, rexdb_tab):
        self.raw_fasta_reads = raw_fasta_reads
        self.outDir = outDir
        self.diamond_path = diamond_path
        self.rexdb_db = self._makeDB(rexdb_fasta) 
        self.rexdb_annotation = rexdb_tab
        self.dic_rexdb_annotation = self.getTE_id_RExDB_dic() #[REXdbIDNumber] = f'{ClassLevel1}_{ClassLevellast}'
        self.outDiamondDir = self._createDic(f'{outDir}/diamond_annotation')
        self.diamond_INFO_per_CLUSTER = {}
        self.reads_in_clusters = {}

    def _createDic(self, directory):
        os.system(f'mkdir {directory}')
        return directory
        
    def _get_target_read_ids(self, tr_seq_reads):
        dic_target = {}
        for seq in SeqIO.parse(tr_seq_reads,'fasta'):
            dic_target[seq.id.split('*')[0]] = 0
        return dic_target
    
    def _makeDB(self,rexdb_fasta):
        os.system(f'diamond makedb  --in {rexdb_fasta} --db {rexdb_fasta}')
        return rexdb_fasta
        
    def getreads_from_cluster(self, cluster, dic_target):
        fout_name = f'{self.outDiamondDir}/clust{cluster}.raw_reads.fasta'
        raw_fasta_reads = self.raw_fasta_reads
        fo_out = open(fout_name, 'w')
        cnt = 0
        for seq in SeqIO.parse(raw_fasta_reads,'fasta'):
             if seq.id in dic_target:
                SeqIO.write(seq, fo_out, 'fasta')
                cnt += 1
        self.reads_in_clusters[f'clust{cluster}'] = cnt
        print(f"Number of reads for cluster clust{cluster} selected - {cnt}")

        fo_out.close()

        return fout_name

    def getTE_id_RExDB_dic(self):
        rexdb_classificatio_dic = {}
        with open(self.rexdb_annotation) as inFile:
            for lines in inFile:
                sp = lines.rstrip().split('\t')
                REXdbIDNumber,  ClassLevel1,  ClassLevellast  = sp[0], sp[1], sp[-1]
                rexdb_classificatio_dic[REXdbIDNumber] = f'{ClassLevel1}_{ClassLevellast}'
        return rexdb_classificatio_dic
    

    def runDiamond(self, cluster):
        tr_seq_reads = f'{self.outDir}/clusters/clust{cluster}.fasta' 
        dic_target = self._get_target_read_ids(tr_seq_reads)
        selected_raw_reads = self.getreads_from_cluster(cluster, dic_target)
        outDiamondDir = self.outDiamondDir
        outTab = f'{outDiamondDir}/{cluster}_vs_rexDB.diamond.tab'
        
        diamon_command = f'diamond blastx --mid-sensitive --long-reads --quiet --db {self.rexdb_db} --query {selected_raw_reads} --range-culling --top 10 --threads 100 --out {outTab} --outfmt 6'
        #print(diamon_command)
        os.system(diamon_command)
        di_hits_per_read_per_classTE, numReadsWithHits, setHits = self.parseDiamonOutTab(outTab)
        self.draw_bar(setHits, cluster)
        self.diamond_INFO_per_CLUSTER[f'clust{cluster}'] = [di_hits_per_read_per_classTE, numReadsWithHits, setHits]
        return [di_hits_per_read_per_classTE, numReadsWithHits, setHits]
    
    def parseDiamonOutTab(self, diamon_blastx_out):
        ''' return list with three items (disctionary (read_id = {'{ClassLevel1_ClassLevellast:[#domain1', domain2],,..}), 
            number of reads with similarity to RExDB and 
            set of unique hits (classification#domain)), 

            see example belowß

            [{'ERR3374012.8976': {'Class_I_Ogre': ['Ty3-INT',
            'Ty3-INT',
            'Ty3-RT',
            'Ty3-RT',
            'Ty3-RT',
            'Ty3-aRH',
            'Ty3-aRH',
            '
            'Ty3-aRH',
            'Ty3-aRH']}},
             1,
             {'Class_I_Ogre#Ty3-INT', 'Class_I_Ogre#Ty3-RT', 'Class_I_Ogre#Ty3-aRH'}]

        '''
        if os.path.getsize(diamon_blastx_out) == 0:
            return [{}, 0, set([])]
        
        TE_id_RExDB_dic = self.dic_rexdb_annotation
        read_diamond_hits = {} # read_id = {'{ClassLevel1_ClassLevellast[#domain1', domain2],,..}
        #return example::

        with open(diamon_blastx_out) as inFile:
            for lines in inFile:
                sp = lines.split('\t')
                read_id, rexdb_id, domain_hit = sp[0], sp[1].split('__')[-1], sp[1].split('__')[0]
                if read_id not in read_diamond_hits:
                    read_diamond_hits[read_id] = {}
                if TE_id_RExDB_dic[rexdb_id].rstrip() not in read_diamond_hits[read_id]:
                    read_diamond_hits[read_id][TE_id_RExDB_dic[rexdb_id]] = []
                read_diamond_hits[read_id][TE_id_RExDB_dic[rexdb_id]].append(domain_hit)
        read_cnt = len(read_diamond_hits)
        hits_list = self._get_domain_list(read_diamond_hits) ## ['dom1:numReads',dom2:numReads']

        return [read_diamond_hits, read_cnt, hits_list]
    
    def _get_domain_list(self,di_hits_per_read_per_classTE):
        rexdb_domain = {}
        for reads in di_hits_per_read_per_classTE:
            per_read_doms = []
            for rexdb_ids in di_hits_per_read_per_classTE[reads]:
                for domains in di_hits_per_read_per_classTE[reads][rexdb_ids]:
                    per_read_doms.append(f'{rexdb_ids}_{domains}')
            for dom_read in set(per_read_doms):
                if dom_read not in rexdb_domain:
                    rexdb_domain[dom_read] = 0
                rexdb_domain[dom_read] += 1
        rexdb_domain = dict(sorted(rexdb_domain.items(), key = lambda x: x[1], reverse = True))
        dom_numReads = [f'{dom}:{rexdb_domain[dom]}' for dom in rexdb_domain]
        return dom_numReads

    def draw_bar(self, setHits, cluster):
        if setHits:
            dom_ids = []
            counts = []
            for doms in setHits:
                dom_id, count = doms.split(":")
                count = int(count)
                dom_ids.append(dom_id)
                counts.append(count)

            plt.style.use('ggplot')
            plt.rcParams["figure.figsize"] = [16, 8]
            font = {'family' : 'sans',
                'weight' : 'bold',
                'size'   : 16}

            plt.rc('font', **font)

            plt.title(f"Reads per TE hit")

            plt.barh(dom_ids, counts)
            # Save the histogram
            plt.savefig(f'{self.outDiamondDir}/clust{cluster}.dom.png')
            plt.clf()