import requests
from ..utils.uniprot import generate_dummy_sets
from ..api.store.client import get_host_param
"""
Run an dummy ORA calculation over a running instance of the pwas service 
"""
def test_pwas_api(pwas_api_host, pwas_api_port, coll_name:str, n_obs:int=1000, f_delta:float=0.05, rh:str="localhost", rp:int=6379, n_top=10):
    uniprot_id_obs, uniprot_id_delta = generate_dummy_sets(coll_name, n_obs, f_delta, rh, rp)

    # Ask Service to compute ora on specified proteome and protein list
    pwas_input = {
        "all_accessions": uniprot_id_obs,
        "taxid" : coll_name,
        "significative_accessions" : uniprot_id_delta,
        "method" : "fisher",
        "maxFreq" : 0.05
        }
    hostname, port = get_host_param()
    ans = requests.post(f"http://{pwas_api_host}:{pwas_api_port}/compute", json = pwas_input)
    if not ans.ok:
        raise ConnectionError("PWAS compute API call failed")
    ora_data = ans.json()
    print(ora_data)
    