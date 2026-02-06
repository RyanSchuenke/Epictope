from epictope.find_executable import find_exe
import os
import subprocess
from Bio.PDB.DSSP import make_dssp_dict
import pandas as pd

def dssp_command(query: str, structure_file: os.PathLike, res_start:int = 1) -> pd.DataFrame:
    """
    Runs dssp on a provided structure file
    
    :param query: accession of the protein dssp is run on
    :type query: str
    :param structure_file: pdb or mmCIF file input to the dssp command
    :type structure_file: os.PathLike
    :param res_start: first residue in the true sequence which appears in the structure file (1 indexed)
    :type res_start: int
    :return: dataframe containing the dssp output
    :rtype: DataFrame
    """
    res_start -= 1
    if not os.path.exists(structure_file):
        raise Exception("Missing structure file for "+query)
    dssp_exe = find_exe("mkdssp")
    out_file = os.path.splitext(structure_file)[0]+".dssp"
    subprocess.run([dssp_exe, structure_file, out_file])
    
    # construct dssp dataframe from saved output file
    dssp = make_dssp_dict(out_file)[0]
    dssp_out = pd.DataFrame(index=range(len(dssp)+res_start), columns=["aa", "structure", "acc", "phi", "psi", "position"])
    for key, value in dssp.items():
        dssp_out.loc[key[1][1]-1+res_start] = list(value[:6])
        dssp_out.at[key[1][1]-1+res_start, "position"]+=res_start
    return dssp_out.set_index(["position", "aa"]).dropna()

def score_ss(dssp:pd.DataFrame, ss_key:dict[int]) -> pd.DataFrame:
    """
    Scores the secondary structure based on the provided ss_key
    
    :param dssp: dataframe containing the mkdssp output
    :type dssp: pd.DataFrame
    :param ss_key: dictionary specifying score for each secondary structure character from 0-1
    :type ss_key: dict[int]
    :return: datafraame with the additional ss_score column
    :rtype: DataFrame
    """
    dssp["ss_score"] = dssp["structure"].map(ss_key)
    return dssp

def rsa(dssp:pd.DataFrame, max_sasa:dict) -> pd.DataFrame:
    """
    Calculates the relative solvent accessibility (RSA)
    
    :param dssp: dataframe containing the mkdssp output
    :type dssp: pd.DataFrame
    :param max_sasa: mapping for the max solvent accessibility for each amino acid
    :type max_sasa: dict
    :return: dataframe with the additional rsa column
    :rtype: DataFrame
    """
    dssp["rsa"] = dssp["acc"] / dssp.index.get_level_values("aa").map(max_sasa)
    return dssp
