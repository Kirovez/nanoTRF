"""
Module to remove unnecessary large files and directories from working directory.
"""
from bin.helpers.help_functions import getLog
import os
import argparse
from Bio import AlignIO
from Bio import SeqIO
import shutil


class Delete_direct():
    def __init__(self, out_dir, dir_clust, dir_canu,dir_reblast, opt_delete,log_file):
        self.out_dir=out_dir
        self.outdir_clust = dir_clust
        self.outdir_canu = dir_canu
        self.outdir_reblast=dir_reblast
        self.del_log = getLog(log_file, "DELETE")


        self.opt_delete = opt_delete
        self.del_dir()
        self.del_log.info("Exit.......\n Finished the work")

    def del_dir(self):
        if self.opt_delete == 'd':
            self.del_log.info("Removing directories has started...")
            #Delete an entire directory tree - ./clust/, ./canu/ and ./ReBlast/
            shutil.rmtree(self.outdir_canu)
            shutil.rmtree(self.outdir_clust)    
            shutil.rmtree(self.outdir_reblast)
            #Delete an TRF html. reports and unnecessary BLAST files
            for file_t in os.listdir(self.out_dir):
                if file_t.endswith('.html') or file_t.endswith('.nhr') or file_t.endswith('.nin') or file_t.endswith('.nsq'):
                    path_t=self.out_dir+file_t
                    os.remove(path_t)
                
            
            
        elif self.opt_delete == 'c':
            self.del_log.info("Directories are not removed")
            pass
        else:
            self.del_log.info("!!!ERROR!!!Parameter does not exist!!!")



