from os.path import exists
from os import PathLike
from yaml import safe_load

DEFAULT_CONFIG = {
    "species" : [
        "bos_taurus", 
        "canis_lupus_familiaris", 
        "gallus_gallus", 
        "homo_sapiens", 
        "mus_musculus", 
        "takifugu_rubripes", 
        "xenopus_tropicalis"
    ],
    "weights" : {
        "h_weight": 1,
        "rsa_weight": 1,
        "ss_weight": 1, 
        "br_weight": 1
    },
    "ss_key" : {
        "G": 0,
        "H": 0,
        "I": 0,
        "E": 0,
        "C": 1,
        "T": 0.5,
        "B": 0.5,
        "S": 0.5,
        "P": 0,
        "-": 1
    },
    "max_sasa" : {
        "A": 121,
        "R": 265,
        "N": 187,
        "D": 187,
        "C": 148,
        "E": 214,
        "Q": 214,
        "G": 97,
        "H": 216,
        "I": 195,
        "L": 191,
        "K": 230,
        "M": 203,
        "F": 228,
        "P": 154,
        "S": 143,
        "T": 163,
        "W": 264,
        "Y": 255,
        "V": 165
    }
}

def load_config(config_path:PathLike = None) -> dict:
    """
    get the config dictionary and update the default config if a custom config is provided
    
    :param config_path: path to custom configuration file
    :type config_path: PathLike
    :return: dictionary containing config parameters
    :rtype: dict
    """
    config = DEFAULT_CONFIG.copy()
    if config_path:
        if exists(config_path):
            print("using custom config file from '"+config_path+"'")
            with open(config_path, 'r') as f:
                config.update(safe_load(f))
        else:
            raise Exception("Config file not found at '"+config_path+"'")
    else: 
        print("using default config")
    return config
