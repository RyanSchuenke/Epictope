from os.path import exists
from os import PathLike
from yaml import safe_load
import logging
logger = logging.getLogger(__name__)

DEFAULT_CONFIG = {
    # define species to run MSA against
    "species" : [
        "bos_taurus", 
        "canis_lupus_familiaris", 
        "gallus_gallus", 
        "homo_sapiens", 
        "mus_musculus", 
        "takifugu_rubripes", 
        "xenopus_tropicalis"
    ],
    # weights for tagging features
    "weights" : {
        "h_weight": 1, # shannon entropy
        "rsa_weight": 1, # solvent accessible surface area
        "ss_weight": 1, # secondary structure
        "br_weight": 1 # disordered binding region
    },
    # value for secondary structures, must be 0-1.
    # each letter refers to a type of secondary structure
    # the number indicates the value or "suitability" for tagging.
    # values should be from 0-1, with higher values indicating greater
    # suitability for tagging.
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
    # reference values for maximum solvent accessibility of amino acids.
    # default values estimate from the following study;
    # https://doi.org/10.1371/journal.pone.0080635
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
            logger.info("using custom config file from '"+config_path+"'")
            with open(config_path, 'r') as f:
                config.update(safe_load(f))
        else:
            logger.error("Config file not found at '"+config_path+"'")
            raise Exception("Config file not found at '"+config_path+"'")
    else: 
        logger.info("using default config")
    return config
