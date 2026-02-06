import os
import subprocess
import tempfile
from epictope.find_executable import find_exe
from Bio import AlignIO, Align

def muscle(query:str, seqs:dict[str], output_folder:os.PathLike) -> Align.MultipleSeqAlignment:
    """
    Infers the multiple sequence alignment of a set of proteins using MUSCLE
    
    :param query: query protein accession
    :type query: str
    :param seqs: dictionary of protein sequences and their accession
    :type seqs: dict[str]
    :param output_folder: path to output folder
    :type output_folder: os.PathLike
    :return: multiple sequence alignment of the inputted proteins
    :rtype: Align.MultipleSeqAlignment
    """
    muscle_exe = find_exe("muscle")
    output_file = os.path.join(output_folder, query)+"_msa.fasta"
    with tempfile.NamedTemporaryFile(mode='w', suffix=".fasta") as tmp:
        for seq in seqs:
            tmp.write(">"+seq+"\n")
            tmp.write(seqs[seq]+"\n")
        tmp.flush()
        
        subprocess.run([muscle_exe, "-align", tmp.name, "-output", output_file])
        with open(output_file, 'r') as file:
            msa = AlignIO.read(file, format="fasta")
    return msa
