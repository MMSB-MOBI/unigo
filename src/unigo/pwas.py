from flask import Flask, request, abort
import requests
from . import vloads as createGOTreeTestUniverseFromAPI
from . import uloads as createGOTreeTestFromAPI

PROXIES = {
  'http': '',
  'https': '',
}


def listen(goApiPort:int):
    global GOPORT
    GOPORT = goApiPort

    app = Flask(__name__)
    app.add_url_rule("/", 'hello', hello)
    app.add_url_rule("/compute", "compute", compute, methods=["POST"])
    return app

def hello():
    return "Hello pwas"

def compute():
    data = request.get_json()
    #Check requested keys
    for key in ["all_accessions", "taxid", "significative_accessions", "method"]:
        if not key in data:
            print(f"ERROR in request : {key} not in posted data. Abort 400")
            abort(400)

    if not data["method"] in ["fisher"]:
        print(f"ERROR : this statistical method ({data['method']}) is not handled. Availables : fisher")
        abort(400)

    all_acc = set(data["all_accessions"])
    significative_acc = set(data["significative_accessions"])
    if not significative_acc.issubset(all_acc):
        print(f"ERROR : significative accessions are not all included in all accessions")
        abort(400)

    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')

    go_url = f"http://127.0.01:{GOPORT}/unigo/{data['taxid']}"
    
    go_resp = requests.get(go_url, proxies = PROXIES)

    if go_resp.status_code != 200:
        print(f"ERROR with go request : {go_url} returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    print(f"GO data loaded from {go_url}")
    tree_universe = createGOTreeTestUniverseFromAPI(go_resp.text)
    n, l, p, p_nr = tree_universe.dimensions
    print(f"Successfully loaded a single {n} nodes and {l} links universal tree")
    print(f"Trying to create Unigo Object using given protein set")
    
    unigoTreeFromAPI = createGOTreeTestFromAPI(go_resp.text, data["all_accessions"])
    x,y = unigoTreeFromAPI.dimensions
    print("Unigo Object successfully buildt w/ following dimensions:")
    print(f"\txpTree => nodes:{x[0]} children_links:{x[1]}, total_protein_occurences:{x[2]}, protein_set:{x[3]}")  
    print(f"\t universeTree => nodes:{y[0]} children_links:{y[1]}, total_protein_occurences:{y[2]}, protein_set:{y[3]}")  

    #Compute fisher stats
    if data["method"] == "fisher":
        print("Computing ORA with fisher")
        rankingsORA = unigoTreeFromAPI.computeORA(data["significative_accessions"], verbose = False)
        
        for i in range(5):

            print(rankingsORA.scoreF[i])
            print(rankingsORA.scoreC[i])

    return "ok"

