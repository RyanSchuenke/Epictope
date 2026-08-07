import requests
import pandas as pd
import os
import sys
import logging
logger = logging.getLogger(__name__)

def config_iupred2a() -> bool:
    """
    Attempts to add the iupred2a directory to path from "IUPRED2A_PATH" environment variable or current working directory
    
    :return: boolean representing if iupred2a package was found or not
    :rtype: bool
    """
    anchor_path = os.getenv("IUPRED2A_PATH")
    
    if anchor_path:
        sys.path.insert(0, anchor_path)
        logger.info("Using anchor2 at '"+anchor_path+"'")
        return True
    else:
        # no environment variable set
        if os.path.exists(os.path.join(os.getcwd(), "iupred2a")):
            sys.path.insert(0, os.path.join(os.getcwd(), "iupred2a"))
            logger.info("Using anchor2 at '"+os.path.join(os.getcwd(), "iupred2a")+"'")
            return True
    # not found in current directory
    return False

def remote_iupred_anchor(uniprot_accession: str) -> pd.DataFrame:
    """
    Retrieve the iupred2/anchor2 data from the iupred2 server
    
    :param uniprot_accession: uniprot accession of the protein to search for
    :type uniprot_accession: str
    :return: dataframe containing iupred2 and anchor2 scores
    :rtype: DataFrame
    """
    iupred_url = "https://iupred2a.elte.hu/iupred2a/anchor/"+uniprot_accession+".json"
    iupred_json = requests.get(iupred_url).json()
    
    anchor_df = pd.DataFrame({"position":range(1,len(iupred_json["sequence"])+1), "aa":list(iupred_json["sequence"]), "iupred2":iupred_json["iupred2"], "anchor2":iupred_json["anchor2"]})
    return anchor_df.set_index(["position", "aa"])

def iupred_anchor(seq: str) -> pd.DataFrame:
    """
    Calculate the iupred2/anchor2 scores with a local installation of iupred2a
    
    :param seq: query protein sequence iupred2/anchor2 calculates on
    :type seq: str
    :return: dataframe containing iupred2 and anchor2 scores
    :rtype: DataFrame
    """
    from iupred2a_lib import iupred, anchor2
    iupred_score = iupred(seq)[0]
    anchor_score = anchor2(seq)
    
    anchor_df = pd.DataFrame({"position":range(1,len(seq)+1),"aa":list(seq),"iupred2a":iupred_score, "anchor2":anchor_score})
    return anchor_df.set_index(["position", "aa"])

def anchor_score(anchor_df:pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the inverted anchor2 score
    
    :param anchor_df: dataframe containing the anchor output data
    :type anchor_df: pd.DataFrame
    :return: anchor dataframe with the additional inv_anchor2 column
    :rtype: DataFrame
    """
    anchor_df["inv_anchor2"] = 1 - anchor_df["anchor2"]
    return anchor_df
