from epictope.setup_folders import setup_folders
from epictope.query_uniprot import query_uniprot
from epictope.fetch_alphafold import fetch_alphafold
from epictope.dssp import dssp_command
from epictope.anchor import iupred_anchor
from epictope.blast import install_db, blast, fetch_seq
from epictope.muscle import muscle
from epictope.shannon import shannon_reshape
from epictope.score import score
from epictope.plot_scores import plot_scores
from epictope.config import load_config
from os import PathLike
from os.path import join
from pandas import DataFrame

def single_score(query:str, config_path:PathLike = None, custom_cif:PathLike = None, res_start:int = 1, graph:bool = False) -> DataFrame:
    if not custom_cif and res_start != 1:
        raise Exception("Cannot set starting residue without a custom structure")
    elif not res_start:
        res_start = 1

    config = load_config(config_path)


    # Setup flders and make blast databases
    folders = setup_folders()

    install_db(species=config["species"], cds_folder=folders["cds_folder"], force=False)

    uniprot_data = query_uniprot(query=query)
    seq = uniprot_data["sequence"]["value"]


    # AlphaFold / DSSP
    if custom_cif:
        print("using custom mmCIF file")
        dssp = dssp_command(query=query, cif_file=custom_cif, res_start=res_start)
        
        if not seq[res_start-1:res_start-1+len(dssp)] == ("".join(dssp.index.get_level_values(1))):
            raise Exception("structure AA sequence does not match protein sequence with starting position "+str(res_start))
    else:
        for cross_ref in uniprot_data["uniProtKBCrossReferences"]:
            if cross_ref["database"] == "AlphaFoldDB":
                alphafold_file = fetch_alphafold(protein_id=cross_ref["id"], model_folder=folders["model_folder"])
                break
        else:
            alphafold_file = fetch_alphafold(query=query, model_folder=folders["model_folder"])
        dssp = dssp_command(query=query, cif_file=alphafold_file)


    # IUPred / Anchor
    anchor_df = iupred_anchor(uniprot_accession=query)


    # Shannon entropy
    blast_hits = {query:seq}
    for species in config["species"]:
        results = blast(seq=seq, db=species, folders=folders)
        blast_hits[species+"_"+results[0][1]] = fetch_seq(seq_id=results[0][1], db=species, cds_folder=folders["cds_folder"])

    alignment = muscle(query=query, seqs=blast_hits, output_folder=folders["output_folder"])

    shannon = shannon_reshape(msa=alignment, query=query, seq_len=len(seq))

    # Calculate final scores
    score_df = score(dssp=dssp, anchor=anchor_df, shannon=shannon, config=config)

    output_file = join(folders["output_folder"], query+"_score")
    score_df.to_csv(output_file+".csv")

    if graph:
        plot_scores(scores_file=output_file+".csv")
    
    return score_df
