from epictope.setup_folders import setup_folders
from epictope.query_uniprot import query_uniprot
from epictope.fetch_alphafold import fetch_alphafold
from epictope.dssp import dssp_command
from epictope.anchor import iupred_anchor
from epictope.blast import install_db, blast, fetch_seq
from epictope.muscle import muscle
from epictope.shannon import shannon_entropy
from epictope.score import score
from epictope.plot_scores import plot_scores
from epictope.config import load_config
from os import PathLike
from os.path import join
from pandas import DataFrame

def single_score(query:str, config_path:PathLike = None, custom_struct:PathLike = None, res_start:int = 1, graph:bool = False) -> DataFrame:
    """
    Function for running the main Epictope pipeline and calculating the 
    "least worst" sites for epitope insertion in a protein sequence.
    
    :param query: Uniprot accession of the protein being scored
    :type query: str
    :param config_path: Path to the config.yml file
    :type config_path: PathLike
    :param custom_struct: Path to a user provided structure file in pdb or cif format for the query protein
    :type custom_struct: PathLike
    :param res_start: 1 indexed starting residue of the protein sequence in the custom_struct file relative to the actual protein sequence
    :type res_start: int
    :param graph: boolean value to determine if the min score should be plotted
    :type graph: bool
    :return: final score dataframe containing shannon entropy, dssp, and iupred2/anchor2 data
    :rtype: DataFrame
    """
    if not custom_struct and res_start != 1:
        raise Exception("Cannot set starting residue without a custom structure")
    elif not res_start:
        res_start = 1

    config = load_config(config_path)


    ## Setup flders and make blast databases
    folders = setup_folders()

    install_db(species=config["species"], cds_folder=folders["cds_folder"], force=False)

    ## Retrieve uniprot data and query sequence
    uniprot_data = query_uniprot(query=query)
    seq = uniprot_data["sequence"]["value"]


    ## AlphaFold / DSSP
    if custom_struct:
        print("using custom structure file")
        dssp = dssp_command(query=query, structure_file=custom_struct, res_start=res_start)
        
        if not seq[res_start-1:res_start-1+len(dssp)] == ("".join(dssp.index.get_level_values(1))):
            raise Exception("structure AA sequence does not match protein sequence with starting position "+str(res_start))
    else:
        # retrieve the alphafold structure by the cross-reference if in uniprot data, otherwise attempt with accession directly
        for cross_ref in uniprot_data["uniProtKBCrossReferences"]:
            if cross_ref["database"] == "AlphaFoldDB":
                alphafold_file = fetch_alphafold(protein_id=cross_ref["id"], model_folder=folders["model_folder"])
                break
        else:
            alphafold_file = fetch_alphafold(query=query, model_folder=folders["model_folder"])
        dssp = dssp_command(query=query, structure_file=alphafold_file)


    ## IUPred / Anchor
    anchor_df = iupred_anchor(uniprot_accession=query)


    ## Shannon entropy
    # 1. find top hit for each of the species
    blast_hits = {query:seq}
    for species in config["species"]:
        results = blast(seq=seq, db=species, folders=folders)
        blast_hits[species+"_"+results[0][1]] = fetch_seq(seq_id=results[0][1], db=species, cds_folder=folders["cds_folder"])

    # 2. create multiple sequence alignment
    alignment = muscle(query=query, seqs=blast_hits, output_folder=folders["output_folder"])

    # 3. calculate shannon entropy
    shannon = shannon_entropy(msa=alignment, query=query, seq_len=len(seq))

    ## Calculate final scores
    score_df = score(dssp=dssp, anchor=anchor_df, shannon=shannon, config=config)

    # Save output file
    output_file = join(folders["output_folder"], query+"_score")
    score_df.to_csv(output_file+".csv")

    # Graph scores
    if graph:
        plot_scores(scores_file=output_file+".csv")
    
    return score_df
