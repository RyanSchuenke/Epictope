import os
import requests
import logging
logger = logging.getLogger(__name__)

def fetch_alphafold(protein_id:str, model_folder:os.PathLike) -> os.PathLike:
    """
    Retrieves the alphafold predicted structure of the query protein from the alphafold database
    
    :param protein_id: protein accession to search for
    :type protein_id: str
    :param model_folder: output folder path for mmCIF file
    :type model_folder: os.PathLike
    :return: path to the downloaded mmCIF file
    :rtype: PathLike
    """
    base_url = "https://alphafold.ebi.ac.uk/files/"
    file_prefix = "AF-"
    version_suffix = "-F1-model_v6.cif"
    
    file_name = "".join([file_prefix, protein_id, version_suffix])
    
    output_path = os.path.join(model_folder, file_name)
    
    if os.path.isfile(output_path):
        logger.info("mmCIF for "+protein_id+" already exists.")
        return output_path
    try:
        r = requests.get(url=base_url+file_name)
        r.raise_for_status()
        with open(output_path, 'wb') as file: 
            file.writelines(r)
    except Exception as err:
        logger.error("Error while downloading: ", protein_id)
        raise err
    return output_path
