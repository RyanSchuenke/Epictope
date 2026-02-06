import requests

def query_uniprot(query:str, fields: list = ["accession", "id", "gene_names", "xref_alphafolddb", "sequence", "organism_name", "organism_id"]) -> dict:
    if not (query or fields):
        raise Exception
    elif len(query) == 0:
        raise Exception("cannot have a 'query' of length zero")
    else: 
        response = requests.get(url = "https://rest.uniprot.org/uniprotkb/"+query+"/", params={"fields":fields})
    return response.json()
