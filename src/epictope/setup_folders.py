import os

def setup_folders() -> list[os.PathLike]:
    folders = {}
    folders["data_folder"] = "data"
    folders["output_folder"] = "outputs"
    folders["model_folder"] = os.path.join(folders["data_folder"],"models")
    folders["cds_folder"] = os.path.join(folders["data_folder"],"CDS")
    folders["temp_folder"] = os.path.join(folders["data_folder"],"temp")
    
    for f in folders:
        os.makedirs(folders[f], exist_ok=True)
    
    return folders
