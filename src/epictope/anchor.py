import requests
import pandas as pd

def iupred_anchor(uniprot_accession: str) -> pd.DataFrame:
    """
    Retrieve the iupred2/anchor2 data from the iupred2 server
    
    :param uniprot_accession: uniprot accession of the protein to search for
    :type uniprot_accession: str
    :return: dataframe containing iupred2 and anchor2 scores
    :rtype: DataFrame
    """
    iupred_url = "https://iupred2a.elte.hu/iupred2a/anchor/"+uniprot_accession+".json"
    iupred_json = requests.get(iupred_url).json()
    
    anchor_df = pd.DataFrame(index=range(len(iupred_json["sequence"])), columns=["position", "aa", "iupred2", "anchor2"])
    for i in range(len(iupred_json["sequence"])):
        anchor_df.loc[i] = [i+1, iupred_json["sequence"][i], iupred_json["iupred2"][i], iupred_json["anchor2"][i]]
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
