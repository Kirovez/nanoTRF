import logging
from pathlib import Path
import os
def checkDir_or_create(dir):

    if not dir.startswith(r"/") and not dir.startswith(r"./"):
        dir = os.path.abspath(os.getcwd()) + "/" + dir

    if not os.path.exists(dir):
        Path(dir).mkdir(parents=True, exist_ok=True)
        logging.info("Working directory was created ....")
    else:
        logging.warning("Working directory existed. It will overwrite output files!!")
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