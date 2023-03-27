import requests
from . import get_host_param

#def unigo_tree_from_api(api_host:str, api_port:int, taxid:int) -> str:
def unigo_tree_from_api(taxid:int) -> str:
    '''Interrogate GO store API and return requests response'''
    hostname, port = get_host_param()
   
    go_url = f"http://{hostname}:{port}/unigos/{taxid}"
    print(f"Fetching GO tree from \"{go_url}\"")
    return requests.get(go_url)

#def unigo_vector_from_api(api_host:str, api_port:int, taxid:int) -> str:
def unigo_vector_from_api(taxid:int) -> str:
    '''Interrogate GO store API and return requests response'''
    hostname, port = get_host_param()
  
    go_url = f"http://{hostname}:{port}/vectors/{taxid}"
    print(f"Fetching GO plain vector from \"{go_url}\"")
    return requests.get(go_url)

#def unigo_culled_from_api(api_host:str, api_port:int, taxid:int, goParameters:{}):
def unigo_culled_from_api(taxid:int, goParameters:{}):
    '''Interrogate GO store API and return requests response'''
    hostname, port = get_host_param()
  
    go_url = f"http://{hostname}:{port}/vectors/{taxid}"
    print(f"Interrogate {go_url} for go culled vector {goParameters}")
    return requests.post(go_url, json=goParameters)

def unigoList(_elem="all"):
    hostname, port = get_host_param()

    d = {}
    for elem in ("trees", "vectors", "culled"):
        if _elem == "all" or elem  == _elem :
            url = f"http://{hostname}:{port}/list/{elem}"
            req = requests.get(url)
            d[elem] = req.json()[elem]
    return d