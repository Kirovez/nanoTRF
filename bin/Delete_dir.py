from bin.helpers.help_functions import getLog
import os
import argparse
from Bio import AlignIO
from Bio import SeqIO
import shutil


class Delete_direct():
    def __init__(self, dir_clust, dir_canu,opt_delete,log_file):
        self.outdir_clust = dir_clust
        self.outdir_canu = dir_canu
        self.del_log = getLog(log_file, "DELETE")


        self.opt_delete = opt_delete
        self.del_dir()
        self.del_log.info("Exit.......\n Finished the work")

    def del_dir(self):
        if self.opt_delete == 'd':
            self.del_log.info("Removing directories has started...")
            shutil.rmtree(self.outdir_canu)
            shutil.rmtree(self.outdir_clust)                                 
            
            
        elif self.opt_delete == 'c':
            self.del_log.info("Directories are not removed")
            pass
        else:
            self.del_log.info("!!!ERROR!!!Parameter does not exist!!!")



