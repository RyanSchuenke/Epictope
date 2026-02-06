from shutil import which
import os

def find_exe(cmd:str="") -> os.PathLike:
    exe = which(cmd)
    if exe:
        return exe
    else: 
        raise Exception("Executable for "+cmd+" not found.")
