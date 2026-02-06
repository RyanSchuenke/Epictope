from epictope.dssp import score_ss, rsa
from epictope.anchor import anchor_score
from epictope.shannon import norm_shannon
import pandas as pd

def score(dssp:pd.DataFrame, anchor:pd.DataFrame, shannon:pd.DataFrame, config:dict) -> pd.DataFrame:
    dssp = score_ss(dssp=dssp, ss_key=config["ss_key"])
    dssp = rsa(dssp)
    anchor = anchor_score(anchor)
    shannon = norm_shannon(shannon)
    
    score_features = ["inv_anchor2","norm_entropy", "rsa", "ss_score"]
    
    score_df = pd.concat([anchor, shannon, dssp], axis=1, join="outer")
    score_df = score_df.infer_objects()
    
    score_df["sum_score"] = score_df[score_features].sum(axis=1, skipna=False)
    score_df["min"] = score_df[score_features].min(axis=1, skipna=False)
    for feature in score_features:
        score_df[f"{feature}_is_min"] = score_df["min"] == score_df[feature]
    
    return score_df
