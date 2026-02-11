from Bio.Align import MultipleSeqAlignment
from math import log2
import pandas as pd

def pos_shannon_entropy(seq:list[str], no_gap:bool = False) -> int:
    """
    Calculates shannon entropy at a single residue
    
    :param seq: List of residues at the same position of a multiple sequence alignment
    :type seq: list[str]
    :param no_gap: Bool for determining if gaps should be included as a unique base for the position
    :type no_gap: bool
    :return: Shannon entropy for the position
    :rtype: int
    """
    if no_gap:
        bases = "".join(seq).replace('-', '')
    else:
        bases = "".join(seq)
    uniq_bases = "".join(set(bases))
    
    n_seqs = len(bases)
    entropy_list = []
    
    for i in uniq_bases:
        # number of times each unique base appears in alignment
        n_i = bases.count(i)
        
        # Pr(base) = # times base appears in alignment / total # sequences in alignment
        prob_i = n_i / n_seqs
        
        # entropy = Pr(base) * log_2(Pr(base))
        entropt_i = prob_i * log2(prob_i)
        entropy_list.append(entropt_i)
    
    # shannon entropy = sum of entopies for each base at the same position of alignment
    res_entropy = -sum(entropy_list)
    return res_entropy if res_entropy else 0

def shannon_entropy(msa:MultipleSeqAlignment, query:str, seq_len:int) -> pd.DataFrame:
    """
    Calculates the shannon entropy for a specific sequence in a multiple sequence alignment
    
    :param msa: Multiple sequence alignments input
    :type msa: MultipleSeqAlignment
    :param query: Query protein of interest in the multiple sequence alignment 
    :type query: str
    :param seq_len: length of query protein for the size of the ouput dataframe
    :type seq_len: int
    :return: DataFrame containing the shannon entropy for each position of the query protein sequence
    :rtype: DataFrame
    """
    for i in range(len(msa)):
        if msa[i].name == query:
            query_index = i
            break
    else:
        raise Exception("Query sequence not found in MSA")
    
    entropy = pd.DataFrame(index=range(seq_len), columns=["position", "aa", "shannon"])
    pos = 0
    
    # Calculate the shannon entropy of each position, only if it is not a gap in the query protein
    for base in range(len(msa[0])):
        if msa[query_index][base] != '-':
            entropy.loc[pos] = (pos+1, msa[query_index][base], pos_shannon_entropy(seq=[seq[base] for seq in msa]))
            pos += 1
    return entropy.set_index(["position", "aa"])

def norm_shannon(shannon: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the normalized shannon entropy for a shannon entropy datafram
    
    :param shannon: DataFrame containing the raw shannon entropy data
    :type shannon: pd.DataFrame
    :return: Input DataFrame with additional normalized shannon entropy column
    :rtype: DataFrame
    """
    # divide by max possible entropy: log_2(20) = 4.321928
    shannon["normalized_entropy"] = shannon["shannon"] / 4.321928
    return shannon
