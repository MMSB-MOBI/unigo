import requests
from ..utils.uniprot import generate_dummy_sets

"""
Run an dummy ORA calculation over a running instance of the pwas service 
"""
def test_pwas_api(api_host:str, api_port:int, coll_name:str, n_obs:int=1000, f_delta:float=0.05, rh:str="localhost", rp:int=6379, n_top=10):
    # 1st try to handshake the pwas server
    try:
        ans = requests.get(f"http://{api_host}:{api_port}/ping")
        if not ans.ok:
            raise ConnectionError("")
    except Exception as e:
        raise ConnectionError("Can't connect to pwas api")
    # Try to Generate uniprot dummy sets
    uniprot_id_obs, uniprot_id_delta = generate_dummy_sets(coll_name, n_obs, f_delta, rh, rp)
    
    # Ask Service to compute ora on specified proteome and protein list
    pwas_input = {
        "all_accessions": uniprot_id_obs,
        "taxid" : coll_name,
        "significative_accessions" : uniprot_id_delta,
        "method" : "fisher",
        "maxFreq" : 0.05
        }
    ans = requests.post(f"http://{api_host}:{api_port}/compute", json = pwas_input)
    if not ans.ok:
        raise ConnectionError("PWAS compute API call failed")
    ora_data = ans.json()
    print(ora_data)
    