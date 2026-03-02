import os
import requests
import logging
logger = logging.getLogger(__name__)

def query_alphafold_uniprot(query:str) -> dict:
    """
    Query Alphafold API for the uniprot query protein sequence and alphafold information
    
    :param query: Uniprot accession of the protein to search for
    :type query: str
    :param fields: list of fields to be included in the response
    :type fields: list
    :return: response dictionary from the uniprot REST api
    :rtype: dict
    """
    if len(query) == 0:
        raise Exception("cannot have a 'query' of length zero")
    else: 
        response = requests.get(url = "https://alphafold.ebi.ac.uk/api/prediction/"+query)
        try:
            response_dict = response.json()[0]
        except:
            logger.error("Could not find "+query+" in alphafold/uniprot database")
            raise Exception("Could not find "+query+" in alphafold/uniprot database")
    return response_dict


def fetch_alphafold(model_url:str, model_folder:os.PathLike) -> os.PathLike:
    """
    Retrieves the alphafold predicted structure of the query protein from the alphafold database
    
    :param protein_id: protein accession to search for
    :type protein_id: str
    :param model_folder: output folder path for mmCIF file
    :type model_folder: os.PathLike
    :return: path to the downloaded mmCIF file
    :rtype: PathLike
    """
    file_name = model_url.split('/')[-1]
    output_path = os.path.join(model_folder, file_name)
    
    if os.path.isfile(output_path):
        logger.info("mmCIF for "+file_name+" already exists.")
        return output_path
    try:
        r = requests.get(url=model_url)
        r.raise_for_status()
        with open(output_path, 'w') as file: 
            file.writelines(r.text)
    except Exception as err:
        logger.error("Error while downloading: ", file_name)
        raise err
    return output_path
