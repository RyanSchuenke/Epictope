import os
import requests

def fetch_alphafold(protein_id:str, model_folder:os.PathLike) -> os.PathLike:
    base_url = "https://alphafold.ebi.ac.uk/files/"
    file_prefix = "AF-"
    version_suffix = "-F1-model_v6.cif"
    
    file_name = "".join([file_prefix, protein_id, version_suffix])
    
    output_path = os.path.join(model_folder, file_name)
    
    if os.path.isfile(output_path):
        print("mmCIF for "+protein_id+" already exists.")
        return output_path
    try:
        r = requests.get(url=base_url+file_name)
        r.raise_for_status()
        with open(output_path, 'wb') as file: 
            file.writelines(r)
    except Exception as err:
        print("Error while downloading: ", protein_id)
        raise err
    return output_path
