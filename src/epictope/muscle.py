import os
import subprocess
import tempfile
from epictope.find_executable import find_exe
from Bio import AlignIO, Align

def muscle(query:str, seqs:dict[str], output_folder:os.PathLike) -> Align.MultipleSeqAlignment:
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
