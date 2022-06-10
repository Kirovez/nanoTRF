import logging
from pathlib import Path
import os



def checkDir_or_create(dir, dir_name = 'Working'):

    if not dir.startswith(r"/") and not dir.startswith(r"./"):
        dir = os.path.abspath(os.getcwd()) + "/" + dir

    if not os.path.exists(dir):
        Path(dir).mkdir(parents=True, exist_ok=True)
        logging.info(f"{dir_name} directory was created ....")
    else:
        logging.warning(f"{dir_name} directory existed. It will overwrite output files!!")
    return dir

def getLog(log_file, log_name):
    LOG = logging.getLogger(log_name)
    LOG.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file, encoding=None, delay=False)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    # add the handlers to the logger
    LOG.addHandler(fh)
    return LOG

def getRExDB():
    rexdb_url = 'http://repeatexplorer.org/repeatexplorer/wp-content/uploads/2018/10/Viridiplantae_v3.0.zip'
    script_directory = os.path.dirname(__file__)
    database_dir = checkDir_or_create(f'{script_directory}/database', dir_name = "RExDB")
    path_to_RExDB_fasta = f'{database_dir}/Viridiplantae_v3.0_ALL_protein-domains.fasta'
    path_to_RExDB_tab = f'{database_dir}/Viridiplantae_v3.0_ALL_classification'
    
    if not os.path.exists(path_to_RExDB_fasta) or not os.path.exists(path_to_RExDB_tab):
        print(f"DOWNLOADING of RExDB files into {database_dir}")
        os.system(f'wget -nc {rexdb_url} -O {database_dir}/Viridiplantae_v3.0_ALL_protein-domains.zip')
        os.system(f'unzip -u {database_dir}/Viridiplantae_v3.0_ALL_protein-domains.zip -d {database_dir}')
        
    assert os.path.exists(path_to_RExDB_fasta), "RExDB fasta file " + path_to_RExDB_fasta + " does not exist!"
    assert os.path.exists(path_to_RExDB_tab), f"RExDB tab file " + path_to_RExDB_tab + " does not exist!"
    
    return(path_to_RExDB_fasta, path_to_RExDB_tab)