import requests
import pandas as pd
import tempfile
import subprocess
import os
import sys
import logging

from epictope.find_executable import find_exe
logger = logging.getLogger(__name__)

def config_phobius() -> bool:
    """
    Attempts to add the phobius directory to path from "PHOBIUS_PATH" environment variable or current working directory
    
    :return: boolean representing if phobius package was found or not
    :rtype: bool
    """
    phobius_path = os.getenv("PHOBIUS_PATH")
    
    if phobius_path:
        sys.path.insert(0, phobius_path)
        logger.info("Using phobius at '"+phobius_path+"'")
        return True
    else:
        # no environment variable set
        if os.path.exists(os.path.join(os.getcwd(), "phobius")):
            os.environ["PATH"] = os.path.join(os.getcwd(), "phobius") + os.pathsep + os.environ.get("PATH", "")
            #sys.path.insert(0, os.path.join(os.getcwd(), "phobius"))
            logger.info("Using phobius at '"+os.path.join(os.getcwd(), "phobius")+"'")
            return True
    # not found in current directory
    return False

def remote_phobius(seq:str) -> pd.DataFrame:
    """
    Predicts transmembrane helices and cellular/extracellular regions usinng the Phobius API
    
    :param seq: query protein sequence
    :type seq: str
    :return: phobius annotation transmembrane prediction by residue
    :rtype: pd.DataFrame
    """
    url = "https://phobius.sbc.su.se/api/predict"
    data = {"sequence":seq}
    response = requests.post(url, json=data).json()
    results = []
    for feature in response['predictions'][0]['features']:
        results.extend(['']*(len(results)+1-feature['start']))
        results.extend([feature['label']]*(feature['stop']-(feature['start']-1)))
    
    phobius_df = pd.DataFrame({"position":range(1,len(seq)+1), "aa":list(seq), "phobius":results})
    phobius_df = phobius_df.set_index(["position", "aa"])
    return phobius_df

def phobius(seq:str) -> pd.DataFrame:
    """
    Predicts transmembrane helices and cellular/extracellular regions usinng a local instalation of Phobius
    
    :param seq: query protein sequence
    :type seq: str
    :return: phobius annotation transmembrane prediction by residue
    :rtype: pd.DataFrame
    """
    try:
        with tempfile.NamedTemporaryFile('w+') as f:
            f.write(">protein\n"+seq)
            f.flush()
            process = subprocess.Popen(
                ["phobius.pl", f.name, '-long'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            output, err = process.communicate()
        output = output.decode('utf8')
        results = []
        for line in output.splitlines()[1:]:
            l = line.split()
            if len(l) > 2:
                if l[1] == "SIGNAL":
                    pass
                elif l[1] == "TRANSMEM":
                    results.extend(['M']*(int(l[3])-(int(l[2])-1)))
                elif l[4].endswith("REGION."):
                    results.extend([l[4][0].lower()]*(int(l[3])-(int(l[2])-1)))
                elif l[4] == "CYTOPLASMIC.":
                    results.extend(['i']*(int(l[3])-(int(l[2])-1)))
                elif l[4] == "NON":
                    results.extend(['o']*(int(l[3])-(int(l[2])-1)))
        phobius_df = pd.DataFrame({"position":range(1,len(seq)+1), "aa":list(seq), "phobius":results})
        phobius_df = phobius_df.set_index(["position", "aa"])
        return phobius_df
    except:
        logger.error(f"Local Phobius failed, switching to api: {err.decode('utf-8')}")
        phobius_df = remote_phobius(seq)
    return phobius_df

def score_phobius(phobius_df:pd.DataFrame) -> pd.DataFrame:
    """
    Scoores residues based on the annotation received from Phobius
    
    :param phobius_df: dataframe containing the phobius output
    :type phobius_df: pd.DataFrame
    :return: input dataframe with additional "phobius_score" column
    :rtype: pd.DataFrame
    """
    feature_key = {
        'i':1, #"inside",
        'o':1, #"outside",
        'O':1, #"long outside",
        'n':0, #"n-terminus",
        'M':0, #"membrane",
        'c':0, #"c-terminus of signal",
        'C':0, #"cleavage site",
        'h':0, #"hydrophobic",
    }
    
    phobius_df["phobius_score"] = phobius_df["phobius"].map(feature_key)
    return phobius_df
