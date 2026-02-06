import requests
import pandas as pd

def iupred_anchor(uniprot_accession: str) -> pd.DataFrame:
    iupred_url = "https://iupred2a.elte.hu/iupred2a/anchor/"+uniprot_accession+".json"
    iupred_json = requests.get(iupred_url).json()
    
    anchor_df = pd.DataFrame(index=range(len(iupred_json["sequence"])), columns=["position", "aa", "iupred2", "anchor2"])
    for i in range(len(iupred_json["sequence"])):
        anchor_df.loc[i] = [i+1, iupred_json["sequence"][i], iupred_json["iupred2"][i], iupred_json["anchor2"][i]]
    return anchor_df.set_index(["position", "aa"])

def anchor_score(anchor_df:pd.DataFrame) -> pd.DataFrame:
    anchor_df["inv_anchor2"] = 1 - anchor_df["anchor2"]
    return anchor_df
