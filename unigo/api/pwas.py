from flask import Flask, abort, jsonify, request
import requests, json
from marshmallow import EXCLUDE
from unigo.api.pwas_compute import compute_over_vector
from .. import uloads as createGOTreeTestFromAPI
from .. import utils
from .io import check_pwas_input_from_route as check_pwas_input
from .data_objects import CulledGoParametersSchema as goParameterValidator
from copy import deepcopy as copy

GOPORT=1234
GOHOST="127.0.0.1"
def listen(goApiHost:str, goApiPort:int, vectorized:bool):
    global GOPORT, GOHOST
    GOPORT = goApiPort
    GOHOST = goApiHost

    app = Flask(__name__)
    app.config['JSON_SORT_KEYS'] = False # To keep dict order in json
    app.add_url_rule("/", 'hello', hello)
    if vectorized:
        print(f"PWAS API vector listening")
        app.add_url_rule("/compute", "computeOverVector", computeOverVector, methods=["POST"])
    else:
        print(f"PWAS API tree listening")
        app.add_url_rule("/compute", "computeOverTree", computeOverTree, methods=["POST"])
    
    app.add_url_rule("/loadVector/<taxid>", loadVector) # WARNING : don't work
    return app

def hello():
    return "Hello pwas"

def computeOverVector():
    data = check_pwas_input()
    return jsonify(compute_over_vector(data, GOHOST, GOPORT, compute_clusters = True))


def computeOverTree():
    data = check_pwas_input() 
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')

    go_resp = utils.unigo_tree_from_api(GOHOST, GOPORT, data["taxid"])

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

def _loadVector(taxid):
    go_resp = utils.unigo_vector_from_api(GOPORT, taxid)
    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    else:
        return go_resp

def loadVector(taxid):
    if _loadVector(taxid):
        return {"ok" : f"Vector loaded for {taxid} taxid"}
    else:
        return {"error" : f"Can't load vector for {taxid} taxid"}




