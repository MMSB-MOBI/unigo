from flask import Flask, request, abort
import requests
from . import vloads as createGOTreeTestUniverseFromAPI
from . import uloads as createGOTreeTestFromAPI
from . import utils


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

    if not utils.check_proteins_subset(data["all_accessions"], data["significative_accessions"]):
        print(f"ERROR : significative accessions are not all included in all accessions")
        abort(400)

    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')

    go_resp = utils.unigo_tree_from_api(GOPORT, data["taxid"])

    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    print("Create Unigo tree")
    unigoTreeFromAPI = createGOTreeTestFromAPI(go_resp.text, data["all_accessions"])
    x,y = unigoTreeFromAPI.dimensions
    print("Unigo Object successfully buildt w/ following dimensions:")
    print(f"\txpTree => nodes:{x[0]} children_links:{x[1]}, total_protein_occurences:{x[2]}, protein_set:{x[3]}")  
    print(f"\t universeTree => nodes:{y[0]} children_links:{y[1]}, total_protein_occurences:{y[2]}, protein_set:{y[3]}")  

    #Compute fisher stats
    if data["method"] == "fisher":
        print("Computing ORA with fisher")
        rankingsORA = unigoTreeFromAPI.computeORA(data["significative_accessions"], verbose = False)
        
        return rankingsORA.json

    return {"not computed": "unavailable stat method"}

