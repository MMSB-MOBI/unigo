from flask import Flask, abort, jsonify
import requests, json

from .. import vloads as createGOTreeTestUniverseFromAPI
from .. import uloads as createGOTreeTestFromAPI
from .. import utils
from .io import checkPwasInput
from ..stat_utils import applyOraToVector, kappaClustering

def listen(goApiPort:int, vectorized:bool):
    global GOPORT
    GOPORT = goApiPort

    app = Flask(__name__)
    app.config['JSON_SORT_KEYS'] = False # To keep dict order in json
    app.add_url_rule("/", 'hello', hello)
    if vectorized:
        print(f"PWAS API vector listening")
        app.add_url_rule("/compute", "computeOverVector", computeOverVector, methods=["POST"])
    else:
        print(f"PWAS API tree listening")
        app.add_url_rule("/compute", "computeOverTree", computeOverTree, methods=["POST"])
    return app

def hello():
    return "Hello pwas"

def computeOverVector():
    data = checkPwasInput() 
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')
    go_resp = utils.unigo_vector_from_api(GOPORT, data["taxid"])

    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    vectorizedProteomeTree = json.loads(go_resp.text)

    res = applyOraToVector(vectorizedProteomeTree, data["all_accessions"], data["significative_accessions"], 0.5)
    Z = {}
    if res:
        Z = kappaClustering(vectorizedProteomeTree["registry"], res)

    complete_results = {
        "list": res,
        "Z" : Z
    }

    print([v["pvalue"] for v in res.values()])

    return jsonify(complete_results)

def computeOverTree():
    data = checkPwasInput() 
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
