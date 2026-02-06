from shutil import which
import os

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
        raise Exception("Executable for "+cmd+" not found.")
