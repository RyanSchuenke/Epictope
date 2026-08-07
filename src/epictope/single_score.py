from epictope.setup_folders import setup_folders
from epictope.fetch_alphafold import query_alphafold_uniprot, fetch_alphafold
from epictope.dssp import dssp_command
from epictope.anchor import remote_iupred_anchor, config_iupred2a, iupred_anchor
from epictope.blast import install_db, blast, fetch_seq
from epictope.muscle import muscle
from epictope.shannon import shannon_entropy
from epictope.score import score
from epictope.plot_scores import plot_scores
from epictope.config import load_config
import os
from Bio import SeqIO
from pandas import DataFrame
import logging
logger = logging.getLogger(__name__)

def single_score(query:str, config_path:os.PathLike = None, custom_struct:os.PathLike = None, plot:bool = False) -> DataFrame:
    """
    Function for running the main Epictope pipeline and calculating the 
    "least worst" sites for epitope insertion in a protein sequence.
    
    :param query: Uniprot accession or fasta file of the protein being scored
    :type query: str
    :param config_path: Path to the config.yml file
    :type config_path: os.PathLike
    :param custom_struct: Path to a user provided structure file in pdb or cif format for the query protein
    :type custom_struct: os.PathLike
    :param plot: boolean value to determine if the min score should be plotted
    :type plot: bool
    :return: final score dataframe containing shannon entropy, dssp, and iupred2/anchor2 data
    :rtype: DataFrame
    """
    config = load_config(config_path)
    
    has_iupred2a = config_iupred2a()


    ## Setup flders and make blast databases
    folders = setup_folders()

    install_db(species=config["species"], cds_folder=folders["cds_folder"], force=False)

    if os.path.exists(query):
        ## Extract query sequence from custom FASTA file 
        if not has_iupred2a:
            logger.error("Cannot use custom fasta sequence file without a local installation of iupred2a")
            raise Exception("Cannot use custom fasta sequence file without a local installation of iupred2a")
        if not custom_struct:
            # maybe allow calling alphafold on the sequence eventually
            logger.error("Cannot use custom fasta sequence file without a custom structure file")
            raise Exception("Cannot use custom fasta sequence file without a custom structure file")
        query_file = query
        seq = str(SeqIO.read(query_file, "fasta").seq)
        query = os.path.splitext(os.path.basename(query))[0]
    else:
        ## Retrieve uniprot data and query sequence
        alphafold_uniprot_data = query_alphafold_uniprot(query=query)
        seq = alphafold_uniprot_data["sequence"]


    ## AlphaFold / DSSP
    if custom_struct:
        logger.info("using custom structure file")
        try:
            struct_seq = str(list(SeqIO.parse(custom_struct, format= "cif-seqres"))[0].seq)
        except:
            try:
                struct_seq = str(list(SeqIO.parse(custom_struct, format= "pdb-atom"))[0].seq)
            except Exception as err:
                logger.error("Failed to open structure file. Please use a .pdb or .cif file")
                raise err
        # Find the starting index of struct_seq in seq
        start_idx = seq.find(struct_seq)
        if start_idx == -1:
            logger.error("The structure sequence is not found in the protein sequence")
            raise Exception("The structure sequence is not found in the protein sequence")
        dssp = dssp_command(structure_file=custom_struct, res_start=start_idx)
    else:
        # retrieve the alphafold structure by the cifUrl
        alphafold_file = fetch_alphafold(model_url = alphafold_uniprot_data["cifUrl"], model_folder=folders["model_folder"])
        dssp = dssp_command(structure_file=alphafold_file)


    ## IUPred / Anchor
    if has_iupred2a:
        # default to local if iupred2a available
        anchor_df = iupred_anchor(seq=seq)
    else:
        # remote lookup via Uniprot accession if no iupred2a provided
        anchor_df = remote_iupred_anchor(uniprot_accession=query)


    ## Shannon entropy
    # 1. find top hit for each of the species
    blast_hits = {query:seq}
    for species in config["species"]:
        results = blast(seq=seq, db=species, folders=folders)
        if results:
            blast_hits[species+"_"+results[0][1]] = fetch_seq(seq_id=results[0][1], db=species, cds_folder=folders["cds_folder"])
        else:
            logger.info("no hits from " + species + " database")

    # 2. create multiple sequence alignment
    alignment = muscle(query=query, seqs=blast_hits, output_folder=folders["output_folder"])

    # 3. calculate shannon entropy
    shannon = shannon_entropy(msa=alignment, query=query, seq_len=len(seq))

    ## Calculate final scores
    score_df = score(dssp=dssp, anchor=anchor_df, shannon=shannon, config=config)

    # Save output file
    output_file = os.path.join(folders["output_folder"], query+"_score")
    score_df.to_csv(output_file+".csv")

    # plot scores
    if plot:
        plot_scores(scores_file=output_file+".csv")
    
    return score_df
