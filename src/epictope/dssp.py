from epictope.find_executable import find_exe
import os
import subprocess
import tempfile
from Bio.PDB.DSSP import make_dssp_dict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
import pandas as pd
import logging
logger = logging.getLogger(__name__)

def dssp_command(structure_file: os.PathLike, res_start:int = 0, save_intermediates: bool = False) -> pd.DataFrame:
    """
    Runs dssp on a provided structure file
    
    :param structure_file: pdb or mmCIF file input to the dssp command
    :type structure_file: os.PathLike
    :param res_start: first residue in the true sequence which appears in the structure file (0 indexed)
    :type res_start: int
    :return: dataframe containing the dssp output
    :rtype: DataFrame
    """
    if not os.path.exists(structure_file):
        logger.error("Missing structure file")
        raise Exception("Missing structure file")
    dssp_exe = find_exe("mkdssp")
    out_file = os.path.splitext(structure_file)[0]+".dssp"
    
    if structure_file.endswith('.cif'):
        # remove checksum from alphafold mmcif files which can cause DSSP error
        with tempfile.NamedTemporaryFile(mode='w+t', suffix=".cif") as structure_in:
            mmcif_dict = MMCIF2Dict(structure_file)
            if "_ma_target_ref_db_details.seq_db_sequence_checksum" in mmcif_dict:
                mmcif_dict.pop("_ma_target_ref_db_details.seq_db_sequence_checksum")
                io=MMCIFIO()
                io.set_dict(mmcif_dict)
                io.save(structure_in.name)
                
            with open(out_file, 'w+t') if save_intermediates else tempfile.NamedTemporaryFile(mode='w+t', suffix=".dssp") as file_out:
                subprocess.run([dssp_exe, structure_in.name, file_out.name])
                dssp = make_dssp_dict(file_out.name)[0]
    else:
        with open(out_file, 'w+t') if save_intermediates else tempfile.NamedTemporaryFile(mode='w+t', suffix=".dssp") as file_out:
            subprocess.run([dssp_exe, structure_file, file_out.name])
            dssp = make_dssp_dict(file_out.name)[0]
    
    # construct dssp dataframe from output file
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
