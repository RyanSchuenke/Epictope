import requests

def query_uniprot(query:str, fields: list = ["accession", "id", "gene_names", "xref_alphafolddb", "sequence", "organism_name", "organism_id"]) -> dict:
    """
    Query uniprot API for the query protein sequence and alphafold information
    
    :param query: Uniprot accession of the protein to search for
    :type query: str
    :param fields: list of fields to be included in the response
    :type fields: list
    :return: response dictionary from the uniprot REST api
    :rtype: dict
    """
    if not (query or fields):
        raise Exception
    elif len(query) == 0:
        raise Exception("cannot have a 'query' of length zero")
    else: 
        response = requests.get(url = "https://rest.uniprot.org/uniprotkb/"+query+"/", params={"fields":fields})
    return response.json()
