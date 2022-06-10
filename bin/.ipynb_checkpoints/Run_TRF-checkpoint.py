import os
import subprocess

def run_TRF(outDir): 
    #running TRF and selecting the sequences from txt.html report
    #running TRF
    consensus_fasta = f'{outDir}/consensus.fasta'
    run_trf=f'trf {consensus_fasta} 2 7 7 80 10 50 2000 -h -d'
    f = open(f'{outDir}/trf.log', 'w') 
    run = subprocess.run(['trf', consensus_fasta, '2','7', '7', '80', '10', '50', '2000', '-h', '-d'],  check=True, stdout=f) #os.system(run_Canu)
    f.close()
    trf_contig1_dic = {} # clustN:[seq_monomer]
    with open(f'{outDir}/trf.parsed','w') as wfile, open('consensus.fasta.2.7.7.80.10.50.2000.dat') as raw_trf:
        seq_id = ''
        for line in raw_trf:
            line = line.rstrip()
            if line.startswith('Sequence'):
                seq_id=line.split(' ')[1]
            if len(line.split(' ')) > 12:
                if seq_id != '':
                    seq_monomer = line.split(' ')[-2]
                    monomer_len = len(seq_monomer)
                    contig_id = seq_id.split('_')[-1]
                    cluster_id = seq_id.split('_')[0]
                    if (contig_id == 'Contig1') or ('artef' in seq_id):
                        if cluster_id not in trf_contig1_dic:
                            trf_contig1_dic[cluster_id] = []
                        trf_contig1_dic[cluster_id].append(seq_monomer)
                    wfile.write(f'{seq_id}\t{seq_monomer}\n')
                    seq_id = ''

    return trf_contig1_dic