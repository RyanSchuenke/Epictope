from Bio.Align import MultipleSeqAlignment
from math import log2
import pandas as pd

def shannon_entropy(seq:list[str], no_gap:bool = False) -> int:
    if no_gap:
        bases = "".join(seq).replace('-', '')
    else:
        bases = "".join(seq)
    uniq_bases = "".join(set(bases))
    
    n_seqs = len(bases)
    entropy_list = []
    
    for i in uniq_bases:
        n_i = bases.count(i)
        prob_i = n_i / n_seqs
        entropt_i = prob_i * log2(prob_i)
        entropy_list.append(entropt_i)
    
    res_entropy = -sum(entropy_list)
    return res_entropy if res_entropy else 0

def shannon_reshape(msa:MultipleSeqAlignment, query:str, seq_len:int) -> pd.DataFrame:
    for i in range(len(msa)):
        if msa[i].name == query:
            query_index = i
            break
    else:
        raise Exception("Query sequence not found in MSA")
    
    entropy = pd.DataFrame(index=range(seq_len), columns=["position", "aa", "shannon"])
    pos = 0
    for base in range(len(msa[0])):
        if msa[query_index][base] != '-':
            entropy.loc[pos] = (pos+1, msa[query_index][base], shannon_entropy(seq=[seq[base] for seq in msa]))
            pos += 1
    return entropy.set_index(["position", "aa"])

def norm_shannon(shannon: pd.DataFrame) -> pd.DataFrame:
    shannon["norm_entropy"] = shannon["shannon"] / 4.321928
    return shannon
