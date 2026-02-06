from epictope.find_executable import find_exe
import os
import subprocess
from Bio.PDB.DSSP import make_dssp_dict
import pandas as pd

def dssp_command(query: str, cif_file: os.PathLike, res_start:int = 1) -> pd.DataFrame:
    res_start -= 1
    if not os.path.exists(cif_file):
        raise Exception("Missing structure file for "+query)
    dssp_exe = find_exe("mkdssp")
    out_file = os.path.splitext(cif_file)[0]+".dssp"
    subprocess.run([dssp_exe, cif_file, out_file])
    
    dssp = make_dssp_dict(out_file)[0]
    dssp_out = pd.DataFrame(index=range(len(dssp)+res_start), columns=["aa", "structure", "acc", "phi", "psi", "position"])
    for key, value in dssp.items():
        dssp_out.loc[key[1][1]-1+res_start] = list(value[:6])
        dssp_out.at[key[1][1]-1+res_start, "position"]+=res_start
    return dssp_out.set_index(["position", "aa"]).dropna()

def score_ss(dssp:pd.DataFrame, ss_key:dict[int]) -> pd.DataFrame:
    dssp["ss_score"] = dssp["structure"].map(ss_key)
    return dssp

def rsa(dssp:pd.DataFrame) -> pd.DataFrame:
    dssp["rsa"] = dssp["acc"] / max(dssp["acc"])
    return dssp
