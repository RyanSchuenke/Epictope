from shutil import which
import os
import logging
logger = logging.getLogger(__name__)

def find_exe(cmd:str="") -> os.PathLike:
    """
    Finds the path to the executable of a program
    
    :param cmd: command to find the executable of
    :type cmd: str
    :return: path to the executable
    :rtype: PathLike
    """
    exe = which(cmd)
    if exe:
        return exe
    else:
        logger.error("Executable for "+cmd+" not found.")
        raise Exception("Executable for "+cmd+" not found.")
