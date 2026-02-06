import os

def setup_folders() -> dict[os.PathLike]:
    """
    Creates the folder structure
    
    :return: dictionary containing the paths of the created folders
    :rtype: dict[PathLike]
    """
    folders = {}
    folders["data_folder"] = "data"
    folders["output_folder"] = "outputs"
    
    # structure and dssp folder
    folders["model_folder"] = os.path.join(folders["data_folder"],"models")
    
    # blast database folder
    folders["cds_folder"] = os.path.join(folders["data_folder"],"CDS") 
    
    for f in folders:
        os.makedirs(folders[f], exist_ok=True)
    
    return folders
