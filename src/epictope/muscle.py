import os
import subprocess
import tempfile
from epictope.find_executable import find_exe
from Bio import AlignIO, Align
import logging
logger = logging.getLogger(__name__)

def muscle(query:str, seqs:dict[str], output_folder:os.PathLike, save_intermediates: bool = False) -> Align.MultipleSeqAlignment:
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
    output_file = os.path.join(output_folder, query)+"_test_2_msa.fasta"

    with tempfile.NamedTemporaryFile(mode='w', suffix=".fasta", delete=False) as tmp_in:
        for seq in seqs:
            tmp_in.write(">"+seq+"\n")
            tmp_in.write(seqs[seq]+"\n")
        tmp_in.flush()
        tmp_name = tmp_in.name
    with open(output_file, 'w+t') if save_intermediates else tempfile.NamedTemporaryFile(mode='w+t', suffix=".fasta") as file_out:
        process = subprocess.Popen(
            [muscle_exe, "-align", tmp_name, "-output", file_out.name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        output, err = process.communicate()
        os.unlink(tmp_name)
        if process.returncode != 0:
            logger.error(f"MUSCLE alignment failed: {err.decode('utf-8')}")
            raise RuntimeError(f"MUSCLE failed: {err.decode('utf-8')}")
        msa = AlignIO.read(file_out.name, format="fasta")
    return msa
